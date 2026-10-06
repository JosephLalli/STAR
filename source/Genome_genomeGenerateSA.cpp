#include "Genome.h"
#include "ErrorWarning.h"
#include <memory>
#define LIBSAIS_OPENMP
#include "libsais/include/libsais64.h"

void Genome::genomeGenerateSA()
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

    if (nGenome>(uint64)numeric_limits<int64_t>::max()/2
            || nGenome>(numeric_limits<size_t>::max()/sizeof(int64_t)-1)/2) {
        exitWithError("EXITING because the libsais suffix array exceeds the addressable size\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
    };
    const uint64 textLength=2*nGenome+1;
    const uint64 packedBytes=!genomeGenerateWriteFiles && P.twoPass.yes ? SApass2.lengthByte : SApass1.lengthByte;
    const uint64 availableRAM=P.limitGenomeGenerateRAM>nG1alloc ? P.limitGenomeGenerateRAM-nG1alloc : 0;
    // The pinned library's per-thread state/cache and aligned buckets fit in this allowance.
    const uint64 scratchBytes=((uint64)P.runThreadN+1)*(1LLU<<19);
    // Sort: 64-bit SA and libsais' worst-case 2*n workspace, with G resident.
    // Pack: the raw and packed arrays coexist after libsais releases its workspace.
    if (scratchBytes>availableRAM || textLength>(availableRAM-scratchBytes)/10
            || packedBytes>availableRAM
            || textLength>(availableRAM-packedBytes)/sizeof(int64_t)) {
        exitWithError("EXITING because --limitGenomeGenerateRAM is too small for the libsais constructor\n"
                      "SOLUTION: increase --limitGenomeGenerateRAM or use --genomeGenerateMethod STAR\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
    };

    unique_ptr<int64_t[]> rawSA;
    try {
        rawSA.reset(new int64_t[textLength]);
    } catch (const bad_alloc &) {
        exitWithError("EXITING because the libsais constructor could not allocate its suffix array\n",
                      std::cerr,P.inOut->logMain,EXIT_CODE_MEMORY_ALLOCATION,P);
    };
    // STAR stops at spacing and orders equal terminated suffixes by position.
    // Distinct boundary bytes in position order preserve that tie rule in libsais.
    for (uint ii=0; ii<nChrReal; ++ii) {
        const uint64 forwardEnd=chrStart[ii]+chrLength[ii];
        const uint64 reverseEnd=2*nGenome-chrStart[nChrReal-1-ii];
        if (G[forwardEnd]!=GENOME_spacingChar || G[reverseEnd]!=GENOME_spacingChar) {
            exitWithError("EXITING because the libsais constructor requires padded chromosome boundaries\n"
                          "SOLUTION: use --genomeGenerateMethod STAR for this genome\n",
                          std::cerr,P.inOut->logMain,EXIT_CODE_INCONSISTENT_DATA,P);
        };
        G[forwardEnd]=GENOME_spacingChar+ii;
        G[reverseEnd]=GENOME_spacingChar+nChrReal+ii;
    };

    const int64_t result=libsais64_omp(reinterpret_cast<const uint8_t*>(G),rawSA.get(),
                                     textLength,0,NULL,P.runThreadN);
    for (uint ii=0; ii<nChrReal; ++ii) {
        G[chrStart[ii]+chrLength[ii]]=GENOME_spacingChar;
        G[2*nGenome-chrStart[ii]]=GENOME_spacingChar;
    };
    if (result!=0) {
        ostringstream errOut;
        errOut << "EXITING because libsais64_omp returned " << result << "\n";
        exitWithError(errOut.str(),std::cerr,P.inOut->logMain,
                      result==-2 ? EXIT_CODE_MEMORY_ALLOCATION : EXIT_CODE_INCONSISTENT_DATA,P);
    };
    if (!genomeGenerateWriteFiles && P.twoPass.yes) {
        SApass2.allocateArray();
        SApass1.pointArray(SApass2.charArray+SApass2.lengthByte-SApass1.lengthByte);
    } else {
        SApass1.allocateArray();
    };
    SA.pointArray(SApass1.charArray+SApass1.lengthByte-SA.lengthByte);
    const uint strandBit=1LLU<<GstrandBit;
    // A/C/G/T-starting suffixes occupy the first nSA ranks: all other symbols exceed 3.
    // Use STAR's serial writer, since neighbouring packed words share bytes.
    for (uint ii=0; ii<nSA; ++ii) {
        const uint position=rawSA[ii];
        SA.writePacked(ii,position<nGenome ? position : ((position-nGenome)|strandBit));
    };
};
