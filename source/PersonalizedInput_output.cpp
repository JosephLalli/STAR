#include "PersonalizedInput.h"
#include "PersonalizedInput_helpers.h"
#include "Genome.h"
#include "Transcriptome.h"
#include "htslib/htslib/bgzf.h"
#include <stdexcept>

using PersonalizedInputHelpers::personalizedOutputPath;

namespace {
const uint64 fastaLineBases=60;
const size_t outputChunkBytes=4ULL<<20;
const size_t gtfQueueBatches=2;

void writeTextFile(const string &path, const string &text)
{
    ofstream out(path+".tmp",ios::binary);
    out << text;
    out.close();
    if (!out || rename((path+".tmp").c_str(),path.c_str())!=0)
        throw runtime_error("cannot write " + path);
}
}

// HTSlib owns BGZF block boundaries, compression and physical .gzi offsets.
class PersonalizedInput::BgzfOutput {
    string path;
    BGZF *file;
    bool indexed;
public:
    uint64 plainBytes=0;
    BgzfOutput(const string &pathIn, bool indexedIn) : path(pathIn),
        file(bgzf_open((path+".tmp").c_str(),"w1")), indexed(indexedIn) {
        if (file==NULL) throw runtime_error("cannot create " + path+".tmp");
        if (indexed && bgzf_index_build_init(file)!=0) {
            bgzf_close(file);
            throw runtime_error("cannot build the BGZF index of " + path);
        };
    }
    ~BgzfOutput() { if (file!=NULL) bgzf_close(file); }
    BgzfOutput(const BgzfOutput &) = delete;
    BgzfOutput &operator=(const BgzfOutput &) = delete;
    void write(const string &text) {
        if (bgzf_write(file,text.data(),text.size())!=static_cast<ssize_t>(text.size()))
            throw runtime_error("cannot write " + path);
        plainBytes+=text.size();
    }
    void finish() {
        if (bgzf_flush(file)!=0 || (indexed && bgzf_index_dump(file,(path+".tmp").c_str(),".gzi")!=0))
            throw runtime_error("cannot finish " + path);
        const int closed=bgzf_close(file);
        file=NULL;
        if (closed!=0 || (indexed && rename((path+".tmp.gzi").c_str(),(path+".gzi").c_str())!=0)
                || rename((path+".tmp").c_str(),path.c_str())!=0)
            throw runtime_error("cannot finish " + path);
    }
};

PersonalizedInput::PersonalizedInput() = default;

PersonalizedInput::~PersonalizedInput()
{
    closeVcfReader();
    closeGtfQueue();
    for (std::thread *worker : {&fastaWriter,&gtfWriter,&transcriptWriter})
        if (worker->joinable()) worker->join();
}

void PersonalizedInput::startOutput(std::thread &worker, OutputState &output,
                                   const std::function<void()> &work)
{
    const auto write=[&output,work]() {
        try { work(); }
        catch (const std::exception &error) { output.error=error.what(); }
    };
    // Keep one core for the producer. Small thread budgets write synchronously.
    if (outputThreads.load()>=requestedThreads-1) {
        write();
    } else {
        ++outputThreads;
        try { worker=std::thread([this,write]() { write(); --outputThreads; }); }
        catch (...) { --outputThreads; throw; }
    };
}

void PersonalizedInput::startFastaOutput(Genome &genome)
{
    if (outputPrefix.empty()) return;
    fastaOutput.path=personalizedOutputPath(outputPrefix,vcfSample,"fa.gz");
    // Native G and chromosome metadata stay valid until the pre-sort output barrier.
    startOutput(fastaWriter,fastaOutput,[this,&genome]() {
        BgzfOutput out(fastaOutput.path,true);
        ostringstream fai;
        string chunk;
        chunk.reserve(outputChunkBytes+fastaLineBases+1);
        for (uint64 contig=0; contig<genome.nChrReal; ++contig) {
            chunk+=">"+genome.chrName[contig]+"\n";
            fai << genome.chrName[contig] << '\t' << genome.chrLength[contig] << '\t'
                << out.plainBytes+chunk.size() << '\t' << fastaLineBases << '\t' << fastaLineBases+1 << '\n';
            const char *bases=genome.G+genome.chrStart[contig];
            for (uint64 line=0; line<genome.chrLength[contig]; line+=fastaLineBases) {
                const uint64 count=min(fastaLineBases,genome.chrLength[contig]-line);
                for (uint64 base=0; base<count; ++base) chunk+="ACGTN"[static_cast<unsigned char>(bases[line+base])];
                chunk+='\n';
                if (chunk.size()>=outputChunkBytes) { out.write(chunk); chunk.clear(); };
            };
        };
        out.write(chunk);
        writeTextFile(fastaOutput.path+".fai",fai.str());
        out.finish();
    });
}

void PersonalizedInput::startTranscriptOutput(Genome &genome)
{
    if (outputPrefix.empty() || !transcriptFasta) return;
    const Transcriptome *annotation=genome.transcriptome;
    if (annotation==NULL) throw runtime_error("the transcript FASTA requires prepared annotation");
    transcriptOutput.path=personalizedOutputPath(outputPrefix,vcfSample,"transcripts.fa.gz");
    startOutput(transcriptWriter,transcriptOutput,[this,&genome,annotation]() {
        BgzfOutput out(transcriptOutput.path,true);
        ostringstream fai;
        string chunk,transcript;
        chunk.reserve(outputChunkBytes+fastaLineBases+1);
        // Native tables give the transcriptome BAM's names, order and exon coordinates.
        for (uint32 index=0; index<annotation->nTr; ++index) {
            transcript.clear();
            const uint64 first=annotation->trExI[index], count=annotation->trExN[index];
            for (uint64 exon=first; exon<first+count; ++exon) {
                const uint64 from=annotation->trS[index]+annotation->exSE[2*exon];
                const uint64 to=annotation->trS[index]+annotation->exSE[2*exon+1];
                for (uint64 base=from; base<=to; ++base) transcript+="ACGTN"[static_cast<unsigned char>(genome.G[base])];
            };
            if (annotation->trStr[index]==2) {
                reverse(transcript.begin(),transcript.end());
                for (char &base : transcript)
                    base=base=='A' ? 'T' : base=='C' ? 'G' : base=='G' ? 'C' : base=='T' ? 'A' : base;
            };
            chunk+=">"+annotation->trID[index]+"\n";
            fai << annotation->trID[index] << '\t' << transcript.size() << '\t'
                << out.plainBytes+chunk.size() << '\t' << fastaLineBases << '\t' << fastaLineBases+1 << '\n';
            for (size_t line=0; line<transcript.size(); line+=fastaLineBases) {
                chunk.append(transcript,line,fastaLineBases);
                chunk+='\n';
                if (chunk.size()>=outputChunkBytes) { out.write(chunk); chunk.clear(); };
            };
        };
        out.write(chunk);
        writeTextFile(transcriptOutput.path+".fai",fai.str());
        out.finish();
    });
}

void PersonalizedInput::queueGtfOutput(const vector<string> &lines)
{
    if (gtfOutput.path.empty()) return;
    string batch;
    for (const string &line : lines) { batch+=line; batch+='\n'; };
    if (!gtfWriter.joinable() && (gtfFile || outputThreads.load()>=requestedThreads-1)) {
        try {
            if (!gtfOutput.error.empty()) return;
            if (!gtfFile) gtfFile.reset(new BgzfOutput(gtfOutput.path,false));
            gtfFile->write(batch);
        } catch (const std::exception &error) { gtfOutput.error=error.what(); }
        return;
    };
    {
        std::unique_lock<std::mutex> lock(gtfQueueMutex);
        gtfQueueReady.wait(lock,[this]() { return gtfQueue.size()<gtfQueueBatches; });
        gtfQueue.push_back(std::move(batch));
    }
    if (!gtfWriter.joinable()) startOutput(gtfWriter,gtfOutput,[this]() { writeGtfOutput(); });
    gtfQueueReady.notify_all();
}

void PersonalizedInput::writeGtfOutput()
{
    unique_ptr<BgzfOutput> out;
    try { out.reset(new BgzfOutput(gtfOutput.path,false)); }
    catch (const std::exception &error) { gtfOutput.error=error.what(); }
    for (;;) {
        string batch;
        {
            std::unique_lock<std::mutex> lock(gtfQueueMutex);
            gtfQueueReady.wait(lock,[this]() { return !gtfQueue.empty() || gtfQueueClosed; });
            if (gtfQueue.empty()) break;
            batch.swap(gtfQueue.front());
            gtfQueue.pop_front();
        }
        gtfQueueReady.notify_all();
        // Drain after failure too, so a bounded producer cannot be left waiting.
        if (out && gtfOutput.error.empty()) {
            try { out->write(batch); }
            catch (const std::exception &error) { gtfOutput.error=error.what(); }
        };
    }
    if (out && gtfOutput.error.empty()) out->finish();
}

void PersonalizedInput::closeGtfQueue()
{
    {
        std::lock_guard<std::mutex> lock(gtfQueueMutex);
        gtfQueueClosed=true;
    }
    gtfQueueReady.notify_all();
}

void PersonalizedInput::finishOutputs()
{
    if (outputsFinished) return;
    closeGtfQueue();
    for (std::thread *worker : {&fastaWriter,&gtfWriter,&transcriptWriter})
        if (worker->joinable()) worker->join();
    if (gtfFile && gtfOutput.error.empty()) {
        try { gtfFile->finish(); }
        catch (const std::exception &error) { gtfOutput.error=error.what(); }
    };
    gtfFile.reset();
    outputsFinished=true;
    for (const OutputState *output : {&fastaOutput,&gtfOutput,&transcriptOutput})
        if (!output->error.empty()) throw runtime_error("personalized output failed: " + output->error);
}
