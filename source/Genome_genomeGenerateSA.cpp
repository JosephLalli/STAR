#include "Genome.h"
#include "ErrorWarning.h"
#include "SjdbClass.h"
#include "sjdbPrepare.h"
#include <memory>
#include <unistd.h>
#define LIBSAIS_OPENMP
#include "libsais/include/libsais64.h"

bool Genome::genomeGenerateSA(SjdbClass &sjdbLoci)
{
    const uint16 byteOrder=1;
    if (pGe.gTypeString!="Full" || pGe.gSAsparseD!=1
            || pGe.gSuffixLengthMax!=numeric_limits<uint>::max() || nChrReal>125
            || *reinterpret_cast<const uint8_t*>(&byteOrder)!=1 || nSA==0) {
        exitWithError("EXITING because --genomeGenerateMethod libsais requires genomeType Full, "
                      "genomeSAsparseD 1, genomeSuffixLengthMax -1, at most 125 contigs, "
                      "little-endian byte order, and A/C/G/T sequence\n"
                      "SOLUTION: use --genomeGenerateMethod STAR for this configuration\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_PARAMETER,P);
    };

    const bool junctionsInSort=P.sjdbInsert.yes && pGe.transform.type==0
            && pGe.sjdbInsertSave!="All" && nChrReal<125;
    const uint64 nGenomeBefore=nGenome;
    const uint64 availableRAM=P.limitGenomeGenerateRAM>nG1alloc ? P.limitGenomeGenerateRAM-nG1alloc : 0;
    uint64 nGsj=0, nSpacers=0, wordWidth=0, digitBase=0, digitCount=0;
    uint64 junctionBufferLength=0;
    unique_ptr<char[]> junctionGsj;

    if (junctionsInSort) {
        if (chrStart[nChrReal]!=nGenomeBefore) {
            exitWithError("EXITING because the libsais junction sort requires the prepared genome to end at the last contig padding\n",
                          std::cerr,P.inOut->logMain,EXIT_CODE_INCONSISTENT_DATA,P);
        };
        P.sjdbInsert.outDir=pGe.gDir;
        P.twoPass.pass2=false;
        junctionBufferLength=2*sjdbLength*sjdbLoci.chr.size()+1;
        if (junctionBufferLength>availableRAM) {
            exitWithError("EXITING because --limitGenomeGenerateRAM is too small for the libsais constructor\n"
                          "SOLUTION: increase --limitGenomeGenerateRAM or use --genomeGenerateMethod STAR\n",
                          std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
        };
        try {
            junctionGsj.reset(new char[junctionBufferLength]);
        } catch (const bad_alloc &) {
            exitWithError("EXITING because the libsais constructor could not allocate its junction sequence buffer\n",
                          std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
        };
        sjdbPrepare(sjdbLoci,P,chrStart[nChrReal],P.sjdbInsert.outDir,*this,junctionGsj.get());
        if (sjdbN>P.limitSjdbInsertNsj) {
            ostringstream errOut;
            errOut << "Fatal LIMIT error: the number of junctions to be inserted on the fly ="<<sjdbN<<" is larger than the limitSjdbInsertNsj="<<P.limitSjdbInsertNsj<<"\n";
            errOut << "SOLUTION: re-run with at least --limitSjdbInsertNsj "<<sjdbN<<"\n";
            exitWithError(errOut.str(),std::cerr,P.inOut->logMain,EXIT_CODE_INPUT_FILES,P);
        };
        nGsj=sjdbN*sjdbLength;
        if (nGsj>0) {
            uint GstrandBit1=(uint)floor(log(nGenomeBefore+nGsj)/log(2))+1;
            if (GstrandBit1<32) GstrandBit1=32;
            if (GstrandBit1>GstrandBit) {
                exitWithError("EXITING because of FATAL ERROR: cannot insert junctions on the fly because of strand GstrandBit problem\n"
                              "SOLUTION: please contact STAR author at https://groups.google.com/forum/#!forum/rna-star\n",
                              std::cerr,P.inOut->logMain,EXIT_CODE_GENOME_FILES,P);
            };
            uint64 junctionBases=0, junctionSpacers=0;
            for (uint64 ii=0; ii<nGsj; ++ii) {
                const uint8 value=static_cast<uint8>(junctionGsj[ii]);
                junctionBases+=value<4;
                junctionSpacers+=value==GENOME_spacingChar;
            };
            nSA+=2*junctionBases;
            nSpacers=2*junctionSpacers+1;
            digitBase=GENOME_spacingChar+1+2*nChrReal;
            digitCount=256-digitBase;
            wordWidth=1;
            for (uint64 rank=nSpacers-1; rank>=digitCount; rank/=digitCount)
                ++wordWidth;
            SA.defineBits(GstrandBit+1,nSA);
        } else {
            junctionGsj.reset();
            junctionBufferLength=0;
        };
    };

    const uint64 nGenomeAfter=nGenomeBefore+nGsj;
    const uint64 tailShift=wordWidth*nSpacers;
    const uint64 signedLimit=numeric_limits<int64_t>::max();
    if (tailShift>=signedLimit || nGenomeAfter>(signedLimit-tailShift-1)/2) {
        exitWithError("EXITING because the libsais suffix array exceeds the addressable size\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
    };
    const uint64 textLength=2*nGenomeAfter+1+tailShift;
    if (textLength>numeric_limits<size_t>::max()/sizeof(int64_t)) {
        exitWithError("EXITING because the libsais suffix array exceeds the addressable size\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
    };

    const uint64 spacerBytes=nSpacers*sizeof(uint64);
    const uint64 packedBytes=!genomeGenerateWriteFiles && P.twoPass.yes ? SApass2.lengthByte : SApass1.lengthByte;
    const uint64 scratchBytes=((uint64)P.runThreadN+1)*(1LLU<<19);
    bool memoryTooSmall=scratchBytes>availableRAM || spacerBytes>availableRAM-scratchBytes;
    if (!memoryTooSmall)
        memoryTooSmall=textLength>(availableRAM-scratchBytes-spacerBytes)/(nGsj>0 ? 11 : 10);
    if (!memoryTooSmall) {
        const uint64 packAvailable=availableRAM-spacerBytes;
        memoryTooSmall=packedBytes>packAvailable
                || textLength>(packAvailable-packedBytes)/sizeof(int64_t);
    };
    if (!memoryTooSmall && nGsj>0)
        memoryTooSmall=junctionBufferLength>availableRAM-spacerBytes
                || textLength>availableRAM-spacerBytes-junctionBufferLength;
    if (memoryTooSmall) {
        exitWithError("EXITING because --limitGenomeGenerateRAM is too small for the libsais constructor\n"
                      "SOLUTION: increase --limitGenomeGenerateRAM or use --genomeGenerateMethod STAR\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
    };

    unique_ptr<uint8_t[]> constructionText;
    vector<uint64> spacerPhysical;
    uint8_t *text=reinterpret_cast<uint8_t*>(G);
    uint64 junctionPhysicalLast=0;
    if (nGsj>0) {
        try {
            constructionText.reset(new uint8_t[textLength]);
            spacerPhysical.reserve(nSpacers);
        } catch (const bad_alloc &) {
            exitWithError("EXITING because the libsais constructor could not allocate its construction text\n",
                          std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
        };
        text=constructionText.get();
        memcpy(text,G,nGenomeBefore);
        uint64 physical=nGenomeBefore;
        const uint8 junctionSpacer=static_cast<uint8>(GENOME_spacingChar+nChrReal);
        auto appendJunctionByte = [&](uint8 value) {
            if (value!=GENOME_spacingChar) {
                text[physical++]=value;
                return;
            };
            const uint64 rank=spacerPhysical.size();
            spacerPhysical.push_back(physical);
            text[physical++]=junctionSpacer;
            uint64 digits=rank;
            for (uint64 digit=wordWidth; digit>0; --digit) {
                text[physical+digit-1]=digitBase+digits%digitCount;
                digits/=digitCount;
            };
            physical+=wordWidth;
        };
        for (uint64 ii=0; ii<nGsj; ++ii)
            appendJunctionByte(junctionGsj[ii]);
        for (uint64 ii=0; ii<nGsj; ++ii) {
            const uint8 value=junctionGsj[nGsj-1-ii];
            appendJunctionByte(value<4 ? 3-value : value);
        };
        appendJunctionByte(GENOME_spacingChar);
        junctionPhysicalLast=spacerPhysical.back()+wordWidth;
        memcpy(text+physical,G+nGenomeBefore+1,nGenomeBefore-1);
        physical+=nGenomeBefore-1;
        text[physical]=GENOME_spacingChar;
        memcpy(G+nGenomeBefore,junctionGsj.get(),nGsj);
        junctionGsj.reset();
    };

    // Rank forward boundaries, junction spacers and reverse boundaries in the
    // order used by STAR when equal suffixes terminate at spacing.
    for (uint ii=0; ii<nChrReal; ++ii) {
        const uint64 forwardEnd=chrStart[ii]+chrLength[ii];
        const uint64 reverseEnd=2*nGenomeAfter-chrStart[nChrReal-1-ii]+tailShift;
        if (text[forwardEnd]!=GENOME_spacingChar || text[reverseEnd]!=GENOME_spacingChar) {
            exitWithError("EXITING because the libsais constructor requires padded chromosome boundaries\n"
                          "SOLUTION: use --genomeGenerateMethod STAR for this genome\n",
                          std::cerr,P.inOut->logMain,EXIT_CODE_INCONSISTENT_DATA,P);
        };
        text[forwardEnd]=GENOME_spacingChar+ii;
        text[reverseEnd]=GENOME_spacingChar+nChrReal+(nGsj>0 ? 1 : 0)+ii;
    };

    const size_t rawSABytes=static_cast<size_t>(textLength)*sizeof(int64_t);
    void *rawSAMapping=mmap(NULL,rawSABytes,PROT_READ|PROT_WRITE,
                            MAP_PRIVATE|MAP_ANON,-1,0);
    if (rawSAMapping==MAP_FAILED) {
        exitWithError("EXITING because the libsais constructor could not allocate its suffix array\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
    };
    auto rawSADeleter = [rawSABytes](int64_t *pointer) {
        if (pointer!=NULL)
            (void)munmap(pointer,rawSABytes);
    };
    unique_ptr<int64_t[],decltype(rawSADeleter)> rawSA(
            static_cast<int64_t*>(rawSAMapping),rawSADeleter);
    const int64_t result=libsais64_omp(text,rawSA.get(),textLength,0,NULL,P.runThreadN);
    if (nGsj==0) {
        for (uint ii=0; ii<nChrReal; ++ii) {
            G[chrStart[ii]+chrLength[ii]]=GENOME_spacingChar;
            G[2*nGenomeBefore-chrStart[ii]]=GENOME_spacingChar;
        };
    };
    if (result!=0) {
        ostringstream errOut;
        errOut << "EXITING because libsais64_omp returned " << result << "\n";
        exitWithError(errOut.str(),std::cerr,P.inOut->logMain,
                      result==-2 ? EXIT_CODE_MEMORY_ALLOCATION : EXIT_CODE_INCONSISTENT_DATA,P);
    };
    if (nGsj>0)
        constructionText.reset();

    if (!genomeGenerateWriteFiles && P.twoPass.yes) {
        SApass2.allocateArray();
        SApass1.pointArray(SApass2.charArray+SApass2.lengthByte-SApass1.lengthByte);
    } else {
        SApass1.allocateArray();
    };
    SA.pointArray(nGsj>0 ? SApass1.charArray : SApass1.charArray+SApass1.lengthByte-SA.lengthByte);
    const uint strandBit=1LLU<<GstrandBit;
    // The first nSA ranks start with A/C/G/T. Convert the physical positions
    // around inserted sentinel words back to STAR's post-insertion coordinates.
    auto packedPosition = [&](uint ii) {
        uint position=rawSA[ii];
        if (nGsj>0 && position>=nGenomeBefore) {
            if (position>junctionPhysicalLast) {
                position-=tailShift;
            } else {
                const uint64 spacersBefore=lower_bound(spacerPhysical.begin(),spacerPhysical.end(),position)-spacerPhysical.begin();
                position-=wordWidth*spacersBefore;
            };
        };
        return position<nGenomeAfter ? position : ((position-nGenomeAfter)|strandBit);
    };
    // Eight packed entries always occupy an integral number of bytes. Each
    // parallel iteration therefore owns a disjoint byte range and can write it
    // without PackedArray::writePacked's unaligned read-modify-write. Bounded
    // blocks let all readers join before complete consumed source pages are released.
    const uint groups=nSA/8;
    const uint groupsPerBlock=(1LLU<<29)/(8*sizeof(int64_t));
    const long pageSizeResult=sysconf(_SC_PAGESIZE);
    const size_t pageBytes=pageSizeResult>0 ? static_cast<size_t>(pageSizeResult) : 0;
    auto releaseRawPages = [&](uint groupBegin, uint groupEnd) {
        if (pageBytes==0)
            return;
        const size_t beginBytes=static_cast<size_t>(groupBegin)*8*sizeof(int64_t);
        const size_t endBytes=static_cast<size_t>(groupEnd)*8*sizeof(int64_t);
        const size_t beginPage=beginBytes/pageBytes+(beginBytes%pageBytes!=0);
        const size_t endPage=endBytes/pageBytes;
        if (endPage>beginPage)
            (void)madvise(reinterpret_cast<char*>(rawSA.get())+beginPage*pageBytes,
                          (endPage-beginPage)*pageBytes,MADV_DONTNEED);
    };
    for (uint groupBegin=0; groupBegin<groups; groupBegin+=groupsPerBlock) {
        const uint groupEnd=groupBegin+min(groupsPerBlock,groups-groupBegin);
        if (P.runThreadN==1) {
            for (uint ii=groupBegin*8; ii<groupEnd*8; ++ii)
                SA.writePacked(ii,packedPosition(ii));
        } else {
            #pragma omp parallel for num_threads(P.runThreadN) schedule(static)
            for (int64 group=(int64)groupBegin; group<(int64)groupEnd; ++group) {
                uint128 buffer=0;
                uint bufferBits=0;
                char *out=SA.charArray+(uint)group*SA.wordLength;
                uint outByte=0;
                for (uint entry=0; entry<8; ++entry) {
                    buffer|=(uint128)packedPosition((uint)group*8+entry)<<bufferBits;
                    bufferBits+=SA.wordLength;
                    while (bufferBits>=8) {
                        out[outByte++]=(char)(buffer&0xff);
                        buffer>>=8;
                        bufferBits-=8;
                    };
                };
            };
        };
        releaseRawPages(groupBegin,groupEnd);
    };
    for (uint ii=groups*8; ii<nSA; ++ii)
        SA.writePacked(ii,packedPosition(ii));

    if (junctionsInSort) {
        if (nGsj>0) {
            nGenome=nGenomeAfter;
            sjGstart=nGenomeBefore;
            SApass1.defineBits(GstrandBit+1,nSA);
            nSAbyte=SA.lengthByte;
        };
        SA.writePacked(nSA,0);
        P.winBinN=nGenome/(1LLU << P.winBinNbits)+1;
    };
    return junctionsInSort;
};
