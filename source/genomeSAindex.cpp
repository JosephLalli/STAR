#include "genomeSAindex.h"
#include "TimeFunctions.h"
#include "SuffixArrayFuns.h"
#include "ErrorWarning.h"
#include <memory>

namespace {

struct IndexBoundary {
    IndexBoundary() : firstIndex(-1), firstIsa(0), lastIndex(0), leadingN(false) {};
    uint firstIndex, firstIsa, lastIndex;
    bool leadingN;
};

void genomeSAindexParallel(char * G, PackedArray & SA, Parameters & P,
                           PackedArray & SAi, Genome &mapGen)
{
    const uint levels=mapGen.pGe.gSAindexNbases;
    const uint noIndex=-1;
    const uint isaStep=mapGen.nSA/(1LLU<<(2*levels))+1;
    const uint threads=min((uint)P.runThreadN,mapGen.nSA);
    const uint span=(mapGen.nSA+threads-1)/threads;
    unique_ptr<uint[]> stage(new uint[mapGen.nSAi]);
    vector<uchar> marked((mapGen.nSAi+7)/8,0);
    vector<IndexBoundary> boundaries(threads*levels);
    int indexError=0;

    #pragma omp parallel for num_threads((int)threads) schedule(static) reduction(|:indexError)
    for (int64 thread=0; thread<(int64)threads; ++thread) {
        const uint first=(uint)thread*span;
        if (first>=mapGen.nSA)
            continue;
        const uint last=min(first+span,mapGen.nSA)-1;
        uint isa=first;
        int iL4;
        uint indFull=funCalcSAiFromSA(G,SA,mapGen,isa,levels,iL4);
        if (thread>0) {
            int previousL4;
            const uint previous=funCalcSAiFromSA(G,SA,mapGen,isa-1,levels,previousL4);
            if (indFull==previous && iL4==previousL4) {
                --isa;
                funSAiFindNextIndex(G,SA,isaStep,isa,indFull,iL4,mapGen);
            };
        };

        vector<uint> lastIndex(levels,noIndex);
        IndexBoundary *boundary=&boundaries[(uint)thread*levels];
        while (isa<=last) {
            for (uint level=0; level<levels; ++level) {
                const uint index=indFull>>(2*(levels-1-level));
                if ((int)level==iL4) {
                    for (uint jj=level; jj<levels; ++jj) {
                        if (lastIndex[jj]==noIndex) {
                            boundary[jj].leadingN=true;
                        } else {
                            const uint cell=mapGen.genomeSAindexStart[jj]+lastIndex[jj];
                            #pragma omp atomic update
                            marked[cell>>3]|=(uchar)(1u<<(cell&7));
                        };
                    };
                    break;
                };

                if (lastIndex[level]==noIndex) {
                    boundary[level].firstIndex=index;
                    boundary[level].firstIsa=isa;
                    lastIndex[level]=index;
                } else if (index>lastIndex[level]) {
                    stage[mapGen.genomeSAindexStart[level]+index]=isa;
                    for (uint jj=lastIndex[level]+1; jj<index; ++jj)
                        stage[mapGen.genomeSAindexStart[level]+jj]=isa|mapGen.SAiMarkAbsentMaskC;
                    lastIndex[level]=index;
                } else if (index<lastIndex[level]) {
                    indexError=1;
                };
            };
            funSAiFindNextIndex(G,SA,isaStep,isa,indFull,iL4,mapGen);
        };
        for (uint level=0; level<levels; ++level) {
            if (lastIndex[level]!=noIndex)
                boundary[level].lastIndex=lastIndex[level];
        };
    };

    vector<uint> previous(levels,noIndex);
    for (uint thread=0; thread<threads; ++thread) {
        for (uint level=0; level<levels; ++level) {
            const IndexBoundary &boundary=boundaries[thread*levels+level];
            if (boundary.leadingN && previous[level]!=noIndex) {
                const uint cell=mapGen.genomeSAindexStart[level]+previous[level];
                marked[cell>>3]|=(uchar)(1u<<(cell&7));
            };
            if (boundary.firstIndex==noIndex)
                continue;
            if (previous[level]==noIndex || boundary.firstIndex>previous[level]) {
                stage[mapGen.genomeSAindexStart[level]+boundary.firstIndex]=boundary.firstIsa;
                const uint firstMissing=previous[level]==noIndex ? 0 : previous[level]+1;
                for (uint jj=firstMissing; jj<boundary.firstIndex; ++jj)
                    stage[mapGen.genomeSAindexStart[level]+jj]=boundary.firstIsa|mapGen.SAiMarkAbsentMaskC;
            } else if (boundary.firstIndex<previous[level]) {
                indexError=1;
            };
            previous[level]=boundary.lastIndex;
        };
    };
    if (indexError!=0) {
        exitWithError("BUG: next index is smaller than previous, EXITING\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_INPUT_FILES,P);
    };

    for (uint level=0; level<levels; ++level) {
        for (uint cell=mapGen.genomeSAindexStart[level]+previous[level]+1;
             cell<mapGen.genomeSAindexStart[level+1]; ++cell)
            stage[cell]=mapGen.nSA|mapGen.SAiMarkAbsentMaskC;
    };

    const uint groups=mapGen.nSAi/8;
    #pragma omp parallel for num_threads((int)threads) schedule(static)
    for (int64 group=0; group<(int64)groups; ++group) {
        uint128 buffer=0;
        uint bufferBits=0;
        char *out=SAi.charArray+(uint)group*SAi.wordLength;
        uint outByte=0;
        for (uint entry=0; entry<8; ++entry) {
            const uint cell=(uint)group*8+entry;
            uint value=stage[cell];
            if ((marked[cell>>3]>>(cell&7))&1u)
                value|=mapGen.SAiMarkNmaskC;
            buffer|=(uint128)value<<bufferBits;
            bufferBits+=SAi.wordLength;
            while (bufferBits>=8) {
                out[outByte++]=(char)(buffer&0xff);
                buffer>>=8;
                bufferBits-=8;
            };
        };
    };
    for (uint cell=groups*8; cell<mapGen.nSAi; ++cell) {
        uint value=stage[cell];
        if ((marked[cell>>3]>>(cell&7))&1u)
            value|=mapGen.SAiMarkNmaskC;
        SAi.writePacked(cell,value);
    };
};

}

void genomeSAindex(char * G, PackedArray & SA, Parameters & P, PackedArray & SAi, Genome &mapGen)
 {
    mapGen.genomeSAindexStart = new uint [mapGen.pGe.gSAindexNbases+1];
    mapGen.genomeSAindexStart[0]=0;
    for (uint ii=1;ii<=mapGen.pGe.gSAindexNbases;ii++) {//L-mer indices starts
        mapGen.genomeSAindexStart[ii] = mapGen.genomeSAindexStart[ii-1] + ( 1LLU<<(2*ii) );
    };
    mapGen.nSAi = mapGen.genomeSAindexStart[mapGen.pGe.gSAindexNbases];

    /* testing
    //     uint* SAi1=new uint[mapGen.nSAi];

    PackedArray SAio;
    SAio.defineBits(mapGen.GstrandBit+3,mapGen.nSAi);
    SAio.allocateArray();
    ifstream oldSAiin("./DirTrue/SAindex");
    oldSAiin.read(SAio.charArray,8*(mapGen.pGe.gSAindexNbases+2));//skip first bytes
    oldSAiin.read(SAio.charArray,SAio.lengthByte);
    oldSAiin.close();
    */


    mapGen.SAiMarkNbit=mapGen.GstrandBit+1;
    mapGen.SAiMarkAbsentBit=mapGen.GstrandBit+2;

    mapGen.SAiMarkNmaskC=1LLU << mapGen.SAiMarkNbit;
    mapGen.SAiMarkNmask=~mapGen.SAiMarkNmaskC;
    mapGen.SAiMarkAbsentMaskC=1LLU << mapGen.SAiMarkAbsentBit;
    mapGen.SAiMarkAbsentMask=~mapGen.SAiMarkAbsentMaskC;


    SAi.defineBits(mapGen.GstrandBit+3,mapGen.nSAi);
    SAi.allocateArray();

    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... generating Suffix Array index\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... generating Suffix Array index\n" <<flush;

    /*testing
    PackedArray SA1=SA;
    uint* ind0=new uint[mapGen.pGe.gSAindexNbases];

    for (uint ii=0; ii<mapGen.pGe.gSAindexNbases; ii++) {
        ind0[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
    };
    uint* SAi1=new uint[mapGen.nSAi];

    for (uint isa=0; isa<mapGen.nSA; isa++) {//for all suffixes
        if (isa%100000000==0) P.inOut->logMain  << isa*100/mapGen.nSA << "% " << flush;

        uint SAstr=SA1[isa];
        bool dirG = (SAstr>>mapGen.GstrandBit) == 0; //forward or reverse strand of the genome
        SAstr &= mapGen.GstrandMask;
        if (!dirG) SAstr=mapGen.nGenome-1-SAstr;

        uint indPref=0;
        for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {//calculate index

            indPref <<= 2;

            uint g1= (uint) G[dirG ? SAstr+iL : SAstr-iL]; //reverese if (-) strand

            if (g1>3) {//if N, this suffix does not belong in SAi
                for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                    SAi1[mapGen.genomeSAindexStart[iL1]+ind0[iL1]] |= mapGen.SAiMarkNmaskC;
                };
                break;
            };

            if (!dirG) g1=3-g1; //complement if (-) strand

            indPref += (uint) g1;

            if ( indPref > ind0[iL] || isa==0 ) {//new && good index, record it
                SAi1[mapGen.genomeSAindexStart[iL]+indPref]=isa;
                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {//index is not present, record to the last present suffix
                    SAi1[mapGen.genomeSAindexStart[iL]+ii] = isa | mapGen.SAiMarkAbsentMaskC;
                };
                ind0[iL]=indPref;
            } else if ( indPref < ind0[iL] ) {
                ostringstream errOut;
                errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            };
        };
    };//for (uint isa=0; isa<mapGen.nSA; isa++)
    */

    uint64 remaining=P.limitGenomeGenerateRAM;
    bool parallelMemory=mapGen.nG1alloc<=remaining;
    if (parallelMemory) remaining-=mapGen.nG1alloc;
    parallelMemory=parallelMemory && mapGen.SApass1.lengthByte<=remaining;
    if (parallelMemory) remaining-=mapGen.SApass1.lengthByte;
    parallelMemory=parallelMemory && SAi.lengthByte<=remaining;
    if (parallelMemory) remaining-=SAi.lengthByte;
    const uint64 temporaryBytes=mapGen.nSAi*sizeof(uint)+(mapGen.nSAi+7)/8;
    parallelMemory=parallelMemory && temporaryBytes<=remaining;

    if (mapGen.pGe.gGenerateMethod=="libsais" && P.runThreadN>1 && parallelMemory) {
        genomeSAindexParallel(G,SA,P,SAi,mapGen);
    } else {
        genomeSAindexChunk(G, SA, P, SAi, 0, SA.length-1, mapGen);
    };

    time(&rawTime);
    P.inOut->logMain    << timeMonthDayTime(rawTime) <<" ... completed Suffix Array index\n" <<flush;
    *P.inOut->logStdOut << timeMonthDayTime(rawTime) <<" ... completed Suffix Array index\n" <<flush;

//     for (uint ii=1;ii<=mapGen.pGe.gSAindexNbases-1;ii++) {//L-mer indices starts
//         cout <<ii<<endl;
//         for (uint jj=mapGen.genomeSAindexStart[ii-1]; jj<mapGen.genomeSAindexStart[ii]; jj++)
//         {
//             if (SAi[jj]!=SAio[jj])
//             {
//                 cout <<ii <<" "<< jj<<" "<<jj-mapGen.genomeSAindexStart[ii-1]<<" "<<SAi[jj]<<" "<<SAio[jj]<<" "<<endl;
//                 sleep(100);
//             };
//         };
//     };


 };


void genomeSAindexChunk(char * G, PackedArray & SA, Parameters & P, PackedArray & SAi, uint iSA1, uint iSA2, Genome &mapGen)
{
    uint* ind0=new uint[mapGen.pGe.gSAindexNbases];

    for (uint ii=0; ii<mapGen.pGe.gSAindexNbases; ii++) {
        ind0[ii]=-1;//this is needed in case "AAA...AAA",i.e. indPref=0 is not present in the genome for some lengths
    };

    PackedArray SAi1;
    SAi1=SAi;
    SAi1.allocateArray();

    uint isaStep=mapGen.nSA/(1llu<<(2*mapGen.pGe.gSAindexNbases))+1;
//     isaStep=8;

    uint isa=iSA1;
    int iL4;
    uint indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4);
    while (isa<=iSA2) {//for all suffixes
        for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {//calculate index

            uint indPref = indFull >> (2*(mapGen.pGe.gSAindexNbases-1-iL));

            if ( (int)iL==iL4 ) {//this suffix contains N and does not belong in SAi
                for (uint iL1=iL; iL1 < mapGen.pGe.gSAindexNbases; iL1++) {
                    SAi.writePacked(mapGen.genomeSAindexStart[iL1]+ind0[iL1],SAi[mapGen.genomeSAindexStart[iL1]+ind0[iL1]] | mapGen.SAiMarkNmaskC);
                };
                break;//break the iL cycle
            };

            if ( indPref > ind0[iL] || isa==0 ) {//new && good index, record it

                SAi.writePacked(mapGen.genomeSAindexStart[iL]+indPref, isa);

                for (uint ii=ind0[iL]+1; ii<indPref; ii++) {//index is not present, record to the last present suffix
                    SAi.writePacked(mapGen.genomeSAindexStart[iL]+ii, isa | mapGen.SAiMarkAbsentMaskC);
                };
                ind0[iL]=indPref;

            } else if ( indPref < ind0[iL] ) {
                ostringstream errOut;
                errOut << "BUG: next index is smaller than previous, EXITING\n" <<flush;
                exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
            };
        };

        //find next index not equal to the current one
        funSAiFindNextIndex(G, SA, isaStep, isa, indFull, iL4, mapGen);//indFull and iL4 have been already defined at the previous step

    };//isa cycle

    for (uint iL=0; iL < mapGen.pGe.gSAindexNbases; iL++) {//fill up unfilled indexes
    	for (uint ii=mapGen.genomeSAindexStart[iL]+ind0[iL]+1; ii<mapGen.genomeSAindexStart[iL+1]; ii++) {
    		SAi.writePacked(ii, mapGen.nSA | mapGen.SAiMarkAbsentMaskC);
    	};
    };

    delete [] ind0;

 };

void funSAiFindNextIndex(char * G, PackedArray & SA, uint isaStep, uint & isa, uint & indFull, int & iL4, Genome &mapGen)
 {
    uint indFullPrev=indFull;
    int iL4prev=iL4;
    isa+=isaStep;
    while (isa<mapGen.nSA && (indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4))==indFullPrev && iL4==iL4prev)
    {//make large step in isa while the indFull/iL4 are still the same
        isa+=isaStep;
    };
    if (isa>=mapGen.nSA)
    {//reached the end of the SA
        indFull=funCalcSAiFromSA(G,SA,mapGen,mapGen.nSA-1,mapGen.pGe.gSAindexNbases,iL4);
        if (indFull==indFullPrev && iL4==iL4prev)
        {
            isa=mapGen.nSA;//no more indices, the last one is equal to the previous
            return;
        };
    };

    {//binary search
        uint i1=isa-isaStep;
        uint i2=min(isa,mapGen.nSA-1);
        while (i1+1<i2)
        {
            isa=i1/2 + i2/2 + (i1%2 + i2%2)/2;
            if ((indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4))==indFullPrev && iL4==iL4prev)
            {
                i1=isa;
            } else
            {
                i2=isa;
            };
        };
        if (isa==i1)
        {
            isa=i2;
            indFull=funCalcSAiFromSA(G,SA,mapGen,isa,mapGen.pGe.gSAindexNbases,iL4);
        };
    };
};
