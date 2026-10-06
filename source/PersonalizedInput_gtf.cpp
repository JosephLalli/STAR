#include "PersonalizedInput.h"

#include "PersonalizedInput_helpers.h"
#include <climits>
#include <stdexcept>

using PersonalizedInputHelpers::splitTabs;

namespace {
const size_t gtfBatchRows = 16384;
const uint64 gtfBatchBytes = 8ULL << 20;
int parseGtfCoordinate(const string &value, const char *label)
{
    if (value.empty()) throw runtime_error(string("malformed GTF ") + label + " coordinate");
    uint64 parsed = 0;
    for (char character : value) {
        if (character < '0' || character > '9')
            throw runtime_error(string("malformed GTF ") + label + " coordinate: " + value);
        if (parsed > (static_cast<uint64>(INT_MAX) - static_cast<uint64>(character - '0')) / 10)
            throw runtime_error(string("GTF ") + label + " coordinate exceeds supported range: " + value);
        parsed = parsed * 10 + static_cast<uint64>(character - '0');
    }
    if (parsed == 0) throw runtime_error(string("GTF ") + label + " coordinate must be positive");
    return static_cast<int>(parsed);
}

string joinTabs(const vector<string> &fields)
{
    string result = fields.at(0);
    for (size_t field = 1; field < fields.size(); ++field)
        result += '\t' + fields[field];
    return result;
}

// The haplotype suffix appended to each gene_id, transcript_id, exon_id and protein_id value, as the
// regex \b(gene_id|transcript_id|exon_id|protein_id) "([^"]+)" -> $1 "$2<suffix>" did: a key counts only
// at the start or after a non-word character, and only with a non-empty quoted value. On all 4,298,549
// matching attributes without a regular-expression pass.
string suffixIdentifiers(const string &attributes, const string &suffix)
{
    static const char *const keys[] = {"gene_id \"", "transcript_id \"", "exon_id \"", "protein_id \""};
    const auto word = [](unsigned char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_';
    };
    string out;
    out.reserve(attributes.size() + 4 * suffix.size());
    size_t copied = 0, i = 0;
    while (i < attributes.size()) {
        size_t close = string::npos;
        if (i == 0 || !word(static_cast<unsigned char>(attributes[i - 1])))
            for (const char *key : keys) {
                const size_t length = strlen(key);
                if (attributes.compare(i, length, key) == 0) {
                    const size_t quote = attributes.find('"', i + length);
                    if (quote != string::npos && quote > i + length) close = quote;
                    break;
                }
            }
        if (close == string::npos) {
            ++i;
            continue;
        }
        out.append(attributes, copied, close - copied);
        out += suffix;
        copied = close;
        i = close + 1;
    }
    out.append(attributes, copied, string::npos);
    return out;
}

} // namespace

PersonalizedInput::GtfResult PersonalizedInput::liftGtfLine(const string &raw) const
{
    GtfResult result;
    try {
        if (raw.empty() || raw[0] == '#') {
            result.lines.push_back(raw);
            return result;
        }
        const vector<string> original = splitTabs(raw);
        if (original.size() != 9 || original[0].empty())
            throw runtime_error("malformed GTF row in " + gtfPath);
        const int start = parseGtfCoordinate(original[3], "start") - 1;
        const int end = parseGtfCoordinate(original[4], "end");
        if (start >= end) throw runtime_error("malformed GTF interval in " + gtfPath);
        result.inputRows = 1;
        result.feature = original[2];
        if (exclude.count(original[0])) {
            result.excludedRows = 1;
            return result;
        }
        if (skipAnnotation.count(original[0])) {
            result.annotationSkippedRows = 1;
            return result;
        }

        const auto source = chainsBySource.find(original[0]);
        if (source == chainsBySource.end())
            throw runtime_error("GTF row has no forward chain for contig: " + original[0]);
        for (const Haplotype &haplotype : source->second) {
            if (start < 0 || end > haplotype.referenceSize)
                throw runtime_error("GTF interval lies outside reference contig: " + original[0]);
            const vector<ChainBlock> &chain = haplotype.chain;
            size_t low = 0, high = chain.size();
            while (low < high) {
                const size_t middle = low + (high - low) / 2;
                if (chain[middle].referenceEnd <= start) low = middle + 1;
                else high = middle;
            }
            const size_t first = low;
            low = 0;
            high = chain.size();
            while (low < high) {
                const size_t middle = low + (high - low) / 2;
                if (chain[middle].referenceStart < end) low = middle + 1;
                else high = middle;
            }
            const size_t afterLast = low;
            if (first >= afterLast || first == chain.size()) {
                ++result.deletedRows;
                continue;
            }

            const ChainBlock &leftBlock = chain[first];
            const ChainBlock &rightBlock = chain[afterLast - 1];
            const int left = max(start, leftBlock.referenceStart);
            const int right = min(end, rightBlock.referenceEnd);
            const int64 queryLeft = static_cast<int64>(leftBlock.queryStart) + left - leftBlock.referenceStart;
            const int64 queryRight = static_cast<int64>(rightBlock.queryStart) + right - rightBlock.referenceStart;
            if (left >= right || queryLeft < 0 || queryLeft >= queryRight || queryRight > haplotype.querySize)
                throw runtime_error("invalid forward-chain projection for GTF row on " + original[0]);

            if (left != start || right != end) ++result.trimmedBoundaryRows;

            vector<string> lifted = original;
            lifted[0] = haplotype.name;
            lifted[3] = to_string(queryLeft + 1);
            lifted[4] = to_string(queryRight);
            if (!haplotype.suffix.empty()) lifted[8] = suffixIdentifiers(lifted[8], haplotype.suffix);
            result.lines.push_back(joinTabs(lifted));
            ++result.outputRows;
        }
    } catch (const std::exception &error) {
        result.error = error.what();
    } catch (...) {
        result.error = "unknown error while lifting personalized GTF row";
    }
    return result;
}

void PersonalizedInput::fillGtfBatch()
{
    if (!gtfReader) {
        if (gtfPath.empty() || gtfPath == "-")
            throw runtime_error("personalized GTF reader requested without a GTF input");
        gtfReader.reset(new PersonalizedInputHelpers::LineReader(gtfPath));
    }

    vector<string> rawLines;
    rawLines.reserve(gtfBatchRows);
    uint64 rawBytes = 0;
    string readError;
    while (rawLines.size() < gtfBatchRows && !gtfInputDone) {
        string raw;
        if (hasPendingGtfRawLine) {
            raw.swap(pendingGtfRawLine);
            hasPendingGtfRawLine = false;
        } else {
            try {
                if (!gtfReader->next(raw)) {
                    gtfInputDone = true;
                    break;
                }
            } catch (const std::exception &error) {
                readError = error.what();
                break;
            } catch (...) {
                readError = "unknown error while reading personalized GTF";
                break;
            }
        }
        const uint64 lineBytes = static_cast<uint64>(raw.size()) + 1;
        if (!rawLines.empty() && lineBytes > gtfBatchBytes - rawBytes) {
            pendingGtfRawLine.swap(raw);
            hasPendingGtfRawLine = true;
            break;
        }
        rawBytes += lineBytes;
        rawLines.push_back(std::move(raw));
        if (rawBytes >= gtfBatchBytes) break;
    }

    if (rawLines.empty()) {
        if (!readError.empty()) pendingGtfError = readError;
        return;
    }
    vector<GtfResult> results(rawLines.size());
    const int budget = availableThreads();
    const int workers = min(budget, static_cast<int>(rawLines.size()));
    if (workers == 1) {
        for (size_t index = 0; index < rawLines.size(); ++index)
            results[index] = liftGtfLine(rawLines[index]);
    } else {
        #pragma omp parallel num_threads(workers) shared(results,rawLines)
        {
            #pragma omp for schedule(static)
            for (int64 index = 0; index < static_cast<int64>(rawLines.size()); ++index)
                results[static_cast<size_t>(index)] = liftGtfLine(rawLines[static_cast<size_t>(index)]);
        }
    }

    for (GtfResult &result : results) {
        if (!result.error.empty()) {
            pendingGtfError = result.error;
            break;
        }
        if (result.inputRows != 0) gtfCounts["input_rows"] += result.inputRows;
        if (result.annotationSkippedRows != 0)
            gtfCounts["annotation_skipped_rows"] += result.annotationSkippedRows;
        if (result.excludedRows != 0) gtfCounts["excluded_rows"] += result.excludedRows;
        if (result.deletedRows != 0) gtfCounts["deleted_rows"] += result.deletedRows;
        if (result.trimmedBoundaryRows != 0)
            gtfCounts["trimmed_boundary_rows"] += result.trimmedBoundaryRows;
        if (result.outputRows != 0) gtfCounts["output_rows"] += result.outputRows;
        if (result.deletedRows != 0)
            gtfCounts["deleted_" + result.feature] += result.deletedRows;
        if (result.outputRows != 0)
            gtfCounts["output_" + result.feature] += result.outputRows;
        for (string &projected : result.lines)
            pendingGtfLines.push_back(std::move(projected));
    }
    queueGtfOutput(pendingGtfLines);
    if (pendingGtfError.empty() && !readError.empty()) pendingGtfError = readError;
}

bool PersonalizedInput::nextGtfLine(string &line)
{
    if (gtfFinished) return false;
    while (true) {
        if (pendingGtfIndex < pendingGtfLines.size()) {
            line = std::move(pendingGtfLines[pendingGtfIndex++]);
            return true;
        }
        pendingGtfLines.clear();
        pendingGtfIndex = 0;
        if (!pendingGtfError.empty()) {
            string error;
            error.swap(pendingGtfError);
            throw runtime_error(error);
        }
        if (gtfInputDone) {
            finishGtf();
            return false;
        }
        fillGtfBatch();
    }
}

void PersonalizedInput::finishGtf()
{
    if (gtfFinished) return;
    gtfFinished = true;
    if (logMain != nullptr) {
        *logMain << "PERSONALIZATION_GTF_COUNTS\n";
        for (const auto &count : gtfCounts)
            *logMain << count.first << "\t" << count.second << "\n";
        *logMain << flush;
    }
    closeGtfQueue();
    gtfReader.reset();
    pendingGtfRawLine.clear();
    hasPendingGtfRawLine = false;
    pendingGtfLines.clear();
    pendingGtfError.clear();
    chainsBySource.clear();
}
