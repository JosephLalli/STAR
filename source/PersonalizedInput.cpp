#include "PersonalizedInput.h"

#include "Genome.h"
#include <cctype>
#include <stdexcept>

#include "PersonalizedInput_helpers.h"

using PersonalizedInputHelpers::personalizedOutputPath;

namespace PersonalizedInputHelpers {
vector<string> splitTabs(const string &line)
{
    vector<string> fields;
    size_t begin = 0;
    while (true) {
        const size_t end = line.find('\t', begin);
        fields.push_back(line.substr(begin, end == string::npos ? string::npos : end - begin));
        if (end == string::npos) return fields;
        begin = end + 1;
    }
}

void mergeCounts(map<string,uint64> &destination, const map<string,uint64> &source)
{
    for (const auto &count : source) destination[count.first] += count.second;
}

string personalizedOutputPath(const string &prefix, const string &sample, const string &kind)
{
    return prefix + sample + ".personalized." + kind;
}
}

namespace {
set<string> parseList(const vector<string> &values)
{
    set<string> result;
    for (const string &value : values)
        if (value != "-") result.insert(value);
    return result;
}

// Whether a donor sample name can name a file: not empty, "." or "..", and without a path separator or whitespace.
bool sampleNamesAFile(const string &sample)
{
    if (sample.empty() || sample == "." || sample == "..") return false;
    for (char c : sample)
        if (c == '/' || isspace(static_cast<unsigned char>(c))) return false;
    return true;
}

} // namespace

std::shared_ptr<PersonalizedInput> PersonalizedInput::build(Genome &genome)
{
    std::shared_ptr<PersonalizedInput> self(new PersonalizedInput);
    self->haploid = parseList(genome.P.personalizationHaploidContigs);
    self->skipAnnotation = parseList(genome.P.personalizationSkipAnnotationContigs);
    self->exclude = parseList(genome.P.personalizationExcludeContigs);
    for (const string &contig : self->exclude)
        if (self->haploid.count(contig) || self->skipAnnotation.count(contig))
            throw runtime_error("contig is both excluded and given another personalization policy: " + contig);
    self->vcfPath = genome.P.personalizationVcf;
    self->vcfSample = genome.P.personalizationSample;
    self->gtfPath = genome.pGe.sjdbGTFfile;
    self->requestedThreads = max(1, genome.P.runThreadN);
    self->logMain = &genome.P.inOut->logMain;
    if (genome.P.personalizationOutputPrefix != "-") self->outputPrefix = genome.P.personalizationOutputPrefix;
    self->transcriptFasta = genome.P.personalizationTranscriptFasta == "Yes";

    if (!self->outputPrefix.empty() && !sampleNamesAFile(self->vcfSample))
        throw runtime_error("the personalization sample cannot name the personalized output files: \"" + self->vcfSample
                            + "\"");

    // Preserve the serial header/index/sample preflight before FASTA input is read.
    self->vcfReader = self->openVcfReader();
    // The GTF writer starts with the first lifted batch, so FASTA preparation keeps every thread.
    if (!self->outputPrefix.empty() && !self->gtfPath.empty() && self->gtfPath != "-")
        self->gtfOutput.path = personalizedOutputPath(self->outputPrefix, self->vcfSample, "gtf.gz");
    return self;
}

int PersonalizedInput::availableThreads() const
{
    return max(1, requestedThreads-outputThreads.load());
}
