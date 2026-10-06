#ifndef STAR_PERSONALIZED_INPUT_H
#define STAR_PERSONALIZED_INPUT_H

#include "IncludeDefine.h"

#include <array>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <iosfwd>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <thread>

class Genome;
namespace PersonalizedInputHelpers { class LineReader; }

// Retains the selected donor-call identity while FASTA is prepared, then the
// forward chain blocks required to project the one-pass GTF reader.
class PersonalizedInput {
public:
    static std::shared_ptr<PersonalizedInput> build(Genome &genome);
    ~PersonalizedInput();

    void appendFasta(Genome &genome);
    bool nextGtfLine(string &line);
    int availableThreads() const;
    void startTranscriptOutput(Genome &genome);
    // Join exporters before STAR mutates or replaces the genome buffer.
    void finishOutputs();

private:
    PersonalizedInput();
    PersonalizedInput(const PersonalizedInput &) = delete;
    PersonalizedInput &operator=(const PersonalizedInput &) = delete;

    struct ChainBlock {
        int referenceStart;
        int referenceEnd;
        int queryStart;
    };
    struct Haplotype {
        string name;
        string suffix;
        int referenceSize;
        int querySize;
        vector<ChainBlock> chain;
    };
    struct SelectedCall {
        int position;
        string reference;
        array<string,2> alternate;
        array<bool,2> applied;
        array<bool,2> sharesAnchor;  // starts on the previous applied call's last REF base
    };
    struct FastaResult {
        vector<string> sequences;
        vector<Haplotype> haplotypes;
        map<string,uint64> counts;
        string error;
    };
    struct GtfResult {
        vector<string> lines;
        string feature;
        string error;
        uint64 inputRows = 0;
        uint64 annotationSkippedRows = 0;
        uint64 excludedRows = 0;
        uint64 deletedRows = 0;
        uint64 trimmedBoundaryRows = 0;
        uint64 outputRows = 0;
    };

    struct OutputState;
    struct VcfReader;
    VcfReader *openVcfReader() const;
    vector<SelectedCall> selectCalls(VcfReader *reader, const string &contig,
                                     const string &reference, map<string,uint64> &counts) const;
    FastaResult prepareContig(const string &source, const string &reference) const;
    GtfResult liftGtfLine(const string &raw) const;
    void fillGtfBatch();
    void closeVcfReader();
    void finishGtf();
    void startFastaOutput(Genome &genome);
    void startOutput(std::thread &worker, OutputState &output, const std::function<void()> &work);
    void queueGtfOutput(const vector<string> &lines);
    void writeGtfOutput();
    void closeGtfQueue();

    struct OutputState {
        string path;
        string error;
    };
    string outputPrefix;
    bool transcriptFasta = false;
    OutputState fastaOutput, gtfOutput, transcriptOutput;
    std::thread fastaWriter, gtfWriter, transcriptWriter;
    class BgzfOutput;
    std::unique_ptr<BgzfOutput> gtfFile;
    std::atomic<int> outputThreads{0};
    std::mutex gtfQueueMutex;
    std::condition_variable gtfQueueReady;
    std::deque<string> gtfQueue;
    bool gtfQueueClosed = false;
    bool outputsFinished = false;

    map<string, vector<Haplotype>> chainsBySource;
    map<string, uint64> callCounts;
    map<string, uint64> gtfCounts;
    set<string> haploid;
    set<string> skipAnnotation;
    set<string> exclude;
    string vcfPath;
    string vcfSample;
    string gtfPath;
    std::unique_ptr<PersonalizedInputHelpers::LineReader> gtfReader;
    string pendingGtfRawLine;
    vector<string> pendingGtfLines;
    size_t pendingGtfIndex = 0;
    string pendingGtfError;
    bool hasPendingGtfRawLine = false;
    bool gtfInputDone = false;
    bool gtfFinished = false;
    int requestedThreads = 1;
    std::ostream *logMain = nullptr;
    VcfReader *vcfReader = nullptr;
};

#endif
