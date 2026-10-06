#include "PersonalizedInput.h"

#include "Genome.h"
#include "SequenceFuns.h"

#include <cctype>
#include <climits>
#include <cstdlib>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <exception>
#include <cstring>
#include <limits>
#include <stdexcept>

#include "PersonalizedInput_helpers.h"

using PersonalizedInputHelpers::splitTabs;
using PersonalizedInputHelpers::mergeCounts;
using PersonalizedInputHelpers::LineReader;

namespace {
const size_t fastaBatchContigs = 12;
const uint64 fastaBatchBytes = 4ULL << 30;
// One contig of an indexed plain FASTA: its .fai name and length, and the byte range of its sequence lines,
// which runs to the next record's header or to the end of the file.
struct IndexedContig {
    string path;
    string name;
    uint64 length = 0;
    uint64 begin = 0;
    uint64 end = 0;
};

struct FileCloser {
    int fd;
    ~FileCloser() { close(fd); }
};

void readExactly(int fd, char *buffer, uint64 bytes, uint64 offset, const string &path)
{
    while (bytes > 0) {
        const ssize_t got = pread(fd, buffer, static_cast<size_t>(min<uint64>(bytes, 1ULL << 30)),
                                  static_cast<off_t>(offset));
        if (got <= 0) throw runtime_error("cannot read FASTA: " + path);
        buffer += got;
        bytes -= static_cast<uint64>(got);
        offset += static_cast<uint64>(got);
    }
}

// Start of the line whose newline is the byte before offset.
uint64 lineStartBefore(int fd, uint64 offset, const string &path)
{
    char block[65536];
    uint64 position = offset - 1;
    while (position > 0) {
        const uint64 begin = position > sizeof(block) ? position - sizeof(block) : 0;
        readExactly(fd, block, position - begin, begin, path);
        for (uint64 i = position - begin; i > 0; --i)
            if (block[i - 1] == '\n') return begin + i;
        position = begin;
    }
    return 0;
}

// The contigs a plain FASTA's .fai names, checked against the file so that a stale or partial index cannot
// change what is read: each indexed header opens its record, the first opens the file, and each record runs
// to the next indexed header. False for gzip, a FIFO, or a missing or empty index: those inputs stream.
bool indexedFastaContigs(const string &path, vector<IndexedContig> &contigs, uint64 &fileBytes)
{
    struct stat info;
    if (stat(path.c_str(), &info) != 0 || !S_ISREG(info.st_mode)) return false;
    ifstream index(path + ".fai");
    if (!index) return false;
    const int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) throw runtime_error("cannot open FASTA: " + path);
    FileCloser closer{fd};
    const uint64 size = static_cast<uint64>(info.st_size);
    unsigned char magic[2] = {0, 0};
    if (size >= 2) readExactly(fd, reinterpret_cast<char *>(magic), 2, 0, path);
    if (magic[0] == 0x1f && magic[1] == 0x8b) return false;

    vector<IndexedContig> found;
    string line;
    uint64 lineNumber = 0;
    while (getline(index, line)) {
        ++lineNumber;
        const vector<string> fields = splitTabs(line);
        const auto number = [&](const string &value) {
            if (value.empty() || value.find_first_not_of("0123456789") != string::npos || value.size() > 19)
                throw runtime_error("malformed FASTA index " + path + ".fai line " + to_string(lineNumber));
            return static_cast<uint64>(stoull(value));
        };
        if (fields.size() < 5 || fields[0].empty())
            throw runtime_error("malformed FASTA index " + path + ".fai line " + to_string(lineNumber));
        IndexedContig contig;
        contig.path = path;
        contig.name = fields[0];
        contig.length = number(fields[1]);
        contig.begin = number(fields[2]);
        found.push_back(contig);
    }
    if (found.empty()) return false;
    for (size_t i = 0; i < found.size(); ++i) {
        IndexedContig &contig = found[i];
        const string mismatch = "FASTA index does not match " + path + " at contig " + contig.name;
        if (contig.begin == 0 || contig.begin > size || (i > 0 && contig.begin <= found[i - 1].begin))
            throw runtime_error(mismatch);
        char newline = 0;
        readExactly(fd, &newline, 1, contig.begin - 1, path);
        if (newline != '\n') throw runtime_error(mismatch);
        const uint64 header = lineStartBefore(fd, contig.begin, path);
        string text(static_cast<size_t>(contig.begin - 1 - header), '\0');
        if (!text.empty()) readExactly(fd, &text[0], text.size(), header, path);
        string name;
        if (text.empty() || text[0] != '>' || !(istringstream(text.substr(1)) >> name) || name != contig.name)
            throw runtime_error(mismatch);
        if (i == 0 ? header != 0 : header < found[i - 1].begin) throw runtime_error(mismatch);
        if (i > 0) found[i - 1].end = header;
    }
    found.back().end = size;
    contigs.insert(contigs.end(), found.begin(), found.end());
    fileBytes += size;
    return true;
}

// One indexed contig's sequence exactly as the streaming reader builds it: bytes below 32 dropped, the rest
// upper-cased. A header line inside the range, or a base count off the index's, means a stale index.
string readIndexedContig(const IndexedContig &contig)
{
    const string mismatch = "FASTA index does not match " + contig.path + " at contig " + contig.name;
    const int fd = open(contig.path.c_str(), O_RDONLY);
    if (fd < 0) throw runtime_error("cannot open FASTA: " + contig.path);
    FileCloser closer{fd};
    string sequence(static_cast<size_t>(contig.end - contig.begin), '\0');
    if (!sequence.empty()) readExactly(fd, &sequence[0], sequence.size(), contig.begin, contig.path);
    size_t kept = 0;
    bool lineStart = true;
    for (size_t i = 0; i < sequence.size(); ++i) {
        const unsigned char base = static_cast<unsigned char>(sequence[i]);
        if (lineStart && base == '>') throw runtime_error(mismatch);
        lineStart = base == '\n';
        if (base >= 32) sequence[kept++] = static_cast<char>(toupper(base));
    }
    if (kept != contig.length) throw runtime_error(mismatch);
    sequence.resize(kept);
    return sequence;
}

} // namespace

void PersonalizedInput::appendFasta(Genome &genome)
{
    uint64 position = 0;
    set<string> seenContigs;
    set<string> excludedSeen;

    const auto pad = [&]() {
        const uint64 bin = genome.genomeChrBinNbases;
        if (bin == 0 || position == UINT64_MAX)
            throw runtime_error("genome coordinate overflow while padding personalized FASTA");
        const uint64 quotient = (position + 1) / bin;
        if (quotient >= UINT64_MAX / bin)
            throw runtime_error("genome coordinate overflow while padding personalized FASTA");
        position = (quotient + 1) * bin;
        if (position > UINT64_MAX - 200)
            throw runtime_error("genome allocation overflow while padding personalized FASTA");
    };
    const auto beginContig = [&](const string &name) {
        if (!genome.chrStart.empty()) {
            genome.chrLength.push_back(position - genome.chrStart.back());
            if (position > 0) pad();
        }
        genome.chrName.push_back(name);
        genome.chrStart.push_back(position);
    };

    struct FastaTask {
        string source;
        string reference;
    };
    struct PreparedContig {
        string source;
        FastaResult result;
        PreparedContig(const string &sourceIn, FastaResult &&resultIn)
            : source(sourceIn), result(std::move(resultIn)) {}
    };
    struct Placement {
        uint64 position;
        const string *sequence;
        Placement(uint64 positionIn, const string *sequenceIn)
            : position(positionIn), sequence(sequenceIn) {}
    };

    vector<FastaTask> tasks;
    vector<PreparedContig> prepared;
    uint64 batchProjectedBytes = 0;
    const uint64 batchLimit = min(fastaBatchBytes, genome.P.limitGenomeGenerateRAM / 16);
    closeVcfReader();

    const auto retainResult = [&](const string &source, FastaResult &&result) {
        if (!result.error.empty()) throw runtime_error(result.error);
        const int outputHaplotypes = haploid.count(source) == 0 ? 2 : 1;
        if (result.sequences.size() != static_cast<size_t>(outputHaplotypes)
                || (!skipAnnotation.count(source)
                    && result.haplotypes.size() != result.sequences.size()))
            throw runtime_error("internal personalized FASTA result mismatch for " + source);
        mergeCounts(callCounts, result.counts);
        prepared.push_back(PreparedContig(source, std::move(result)));
    };

    const auto flushBatch = [&]() {
        if (tasks.empty()) return;
        vector<FastaResult> results(tasks.size());
        const int workers = min(requestedThreads, static_cast<int>(tasks.size()));
        if (workers == 1) {
            for (size_t index = 0; index < tasks.size(); ++index)
                results[index] = prepareContig(tasks[index].source, tasks[index].reference);
        } else {
            #pragma omp parallel num_threads(workers) shared(results,tasks)
            {
                #pragma omp for schedule(static)
                for (int64 index = 0; index < static_cast<int64>(tasks.size()); ++index)
                    results[static_cast<size_t>(index)] = prepareContig(
                        tasks[static_cast<size_t>(index)].source,
                        tasks[static_cast<size_t>(index)].reference);
            }
        }

        for (size_t index = 0; index < tasks.size(); ++index)
            retainResult(tasks[index].source, std::move(results[index]));
        tasks.clear();
        batchProjectedBytes = 0;
    };

    const auto queueContig = [&](const string &source, string &reference) {
        if (exclude.count(source)) {
            excludedSeen.insert(source);
            reference.clear();
            return;
        }
        if (reference.empty() || reference.size() > static_cast<size_t>(INT_MAX)) {
            flushBatch();
            throw runtime_error("reference contig must have a positive signed-32-bit length: " + source);
        }
        if (seenContigs.count(source)) {
            flushBatch();
            throw runtime_error("duplicate FASTA contig name in personalized input: " + source);
        }
        const uint64 copies = haploid.count(source) == 0 ? 3 : 2;
        const uint64 projected = copies * static_cast<uint64>(reference.size());
        const bool exceedsLimit = projected > batchLimit
            || batchProjectedBytes > batchLimit - projected;
        if (!tasks.empty() && (tasks.size() >= fastaBatchContigs || exceedsLimit)) flushBatch();
        seenContigs.insert(source);
        tasks.push_back({source, std::move(reference)});
        batchProjectedBytes += projected;
        reference.clear();
    };

    // A plain FASTA with a .fai is read by contig from every worker at once: each reads and personalizes its
    // contigs, largest first, and the genome is laid out once they finish.
    vector<IndexedContig> indexed;
    uint64 indexedFileBytes = 0;
    bool indexedIngest = !genome.pGe.gFastaFiles.empty();
    for (const string &path : genome.pGe.gFastaFiles)
        indexedIngest = indexedIngest && indexedFastaContigs(path, indexed, indexedFileBytes);
    if (indexedIngest) {
        vector<IndexedContig> queued;
        string queueError;
        for (const IndexedContig &contig : indexed) {
            if (exclude.count(contig.name)) {
                excludedSeen.insert(contig.name);
                continue;
            }
            if (contig.length == 0 || contig.length > static_cast<uint64>(INT_MAX)) {
                queueError = "reference contig must have a positive signed-32-bit length: " + contig.name;
                break;
            }
            if (seenContigs.count(contig.name)) {
                queueError = "duplicate FASTA contig name in personalized input: " + contig.name;
                break;
            }
            seenContigs.insert(contig.name);
            queued.push_back(contig);
        }
        vector<size_t> order(queued.size());
        for (size_t index = 0; index < order.size(); ++index) order[index] = index;
        stable_sort(order.begin(), order.end(),
                    [&](size_t a, size_t b) { return queued[a].length > queued[b].length; });
        vector<FastaResult> results(queued.size());
        const int workers = max(1, min(requestedThreads, static_cast<int>(queued.size())));
        #pragma omp parallel num_threads(workers)
        {
            #pragma omp for schedule(dynamic,1)
            for (int64 slot = 0; slot < static_cast<int64>(order.size()); ++slot) {
                const size_t index = order[static_cast<size_t>(slot)];
                try {
                    results[index] = prepareContig(queued[index].name, readIndexedContig(queued[index]));
                } catch (const std::exception &error) {
                    results[index].error = error.what();
                }
            }
        }
        for (size_t index = 0; index < queued.size(); ++index)
            retainResult(queued[index].name, std::move(results[index]));
        if (!queueError.empty()) throw runtime_error(queueError);
    }

    for (const string &path : genome.pGe.gFastaFiles) {
        if (indexedIngest) break;
        LineReader input(path);
        string line;
        string source;
        string reference;
        bool firstLine = true;
        while (true) {
            bool haveLine = false;
            try {
                haveLine = input.next(line);
            } catch (...) {
                flushBatch();
                throw;
            }
            if (!haveLine) break;
            if (firstLine && (line.empty() || line[0] != '>')) {
                flushBatch();
                throw runtime_error("FASTA must start with a header: " + path);
            }
            firstLine = false;
            if (!line.empty() && line[0] == '>') {
                if (!source.empty()) queueContig(source, reference);
                istringstream header(line.substr(1));
                if (!(header >> source)) {
                    flushBatch();
                    throw runtime_error("empty FASTA contig name: " + path);
                }
                reference.clear();
                continue;
            }
            if (source.empty()) {
                flushBatch();
                throw runtime_error("FASTA sequence appears before a header: " + path);
            }
            for (char base : line)
                if (static_cast<unsigned char>(base) >= 32)
                    reference.push_back(static_cast<char>(toupper(static_cast<unsigned char>(base))));
        }
        if (firstLine) {
            flushBatch();
            throw runtime_error("empty FASTA input: " + path);
        }
        if (source.empty()) {
            flushBatch();
            throw runtime_error("FASTA contains no contig header: " + path);
        }
        queueContig(source, reference);
        flushBatch();
    }

    // An excluded name the FASTA never held is a mistaken layout, not a no-op.
    for (const string &contig : exclude)
        if (!excludedSeen.count(contig))
            throw runtime_error("--personalizationExcludeContigs names a contig absent from the FASTA: " + contig);

    vector<Placement> placements;
    for (PreparedContig &contig : prepared) {
        const int outputHaplotypes = haploid.count(contig.source) == 0 ? 2 : 1;
        for (int haplotype = 0; haplotype < outputHaplotypes; ++haplotype) {
            const string suffix = outputHaplotypes == 1
                ? string() : (haplotype == 0 ? "_L" : "_R");
            beginContig(contig.source + suffix);
            const string &sequence = contig.result.sequences[static_cast<size_t>(haplotype)];
            if (sequence.size() > static_cast<size_t>(INT_MAX)
                    || position > UINT64_MAX - static_cast<uint64>(sequence.size()) - 200)
                throw runtime_error("personalized FASTA sequence/allocation exceeds supported range");
            placements.push_back(Placement(position, &sequence));
            position += sequence.size();
        }
        if (!contig.result.haplotypes.empty())
            chainsBySource[contig.source] = std::move(contig.result.haplotypes);
    }

    if (genome.chrStart.empty() || position == 0)
        throw runtime_error("FASTA contains no sequence");
    genome.chrLength.push_back(position - genome.chrStart.back());
    pad();
    if (position > numeric_limits<uint64>::max() / 2 - 100
            || position > static_cast<uint64>(numeric_limits<size_t>::max() / 2 - 100))
        throw runtime_error("genome allocation overflow for the personalized FASTA");
    const uint64 nativeBytes = 2 * (position + 100);
    if (nativeBytes > numeric_limits<uint64>::max() - nativeBytes / 3)
        throw runtime_error("genome allocation overflow for the personalized FASTA");

    genome.nGenome = position;
    genome.nChrReal = genome.chrName.size();
    genome.chrStart.push_back(position);
    for (uint64 chromosome = 0; chromosome < genome.nChrReal; ++chromosome)
        genome.chrNameIndex[genome.chrName[chromosome]] = chromosome;
    genome.genomeSequenceAllocate(genome.nGenome, genome.nG1alloc, genome.G, genome.G1);

    vector<unsigned char> copyFailed(placements.size(), 0);
    #pragma omp parallel for num_threads(max(1, requestedThreads)) schedule(dynamic,1)
    for (int64 piece = 0; piece < static_cast<int64>(placements.size()); ++piece) {
        const Placement &placement = placements[static_cast<size_t>(piece)];
        const string &sequence = *placement.sequence;
        const uint written = convertNucleotidesToNumbersRemoveControls(
            sequence.data(), genome.G + placement.position, static_cast<uint>(sequence.size()));
        if (written != sequence.size()) copyFailed[static_cast<size_t>(piece)] = 1;
    }
    for (unsigned char failed : copyFailed)
        if (failed != 0)
            throw runtime_error("unexpected control character after personalized FASTA normalization");

    vector<Placement>().swap(placements);
    vector<PreparedContig>().swap(prepared);

    for (const string &contig : exclude)
        genome.P.inOut->logMain << "STAR_PERSONALIZATION excluded_contig=" << contig << "\n";
    genome.P.inOut->logMain << "PERSONALIZATION_CALL_COUNTS\n";
    for (const auto &count : callCounts)
        genome.P.inOut->logMain << count.first << "\t" << count.second << "\n";
    genome.P.inOut->logMain << flush;
    startFastaOutput(genome);
}
