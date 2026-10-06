#include "GTF.h"
#include "PersonalizedInput.h"

#include "ErrorWarning.h"
#include "streamFuns.h"
#include "TimeFunctions.h"

GTF::GTF(Genome &genome, Parameters &P, const string &dirOut, SjdbClass &sjdbLoci) 
                 : genome(genome), P(P), dirOut(dirOut), sjdbLoci(sjdbLoci), superTrome(P)
 {//initialize; load gtf file; returns number of added junctions
                     
    if (genome.sjdbOverhang==0 || genome.pGe.sjdbGTFfile=="-") {//no GTF
        gtfYes=false;
        return;
    };
    gtfYes=true;
    mkdir(dirOut.c_str(), P.runDirPerm);
    
    time_t rawTime;
    time(&rawTime);
    P.inOut->logMain     << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;
    *P.inOut->logStdOut  << timeMonthDayTime(rawTime) <<" ..... processing annotations GTF\n" <<flush;

    std::map <string,uint64> transcriptIDnumber, geneIDnumber;    
    
    ifstream sjdbStreamIn;
    if (!genome.personalizedInput) sjdbStreamIn.open(genome.pGe.sjdbGTFfile.c_str());
    if (!genome.personalizedInput && sjdbStreamIn.fail()) {
        ostringstream errOut;
        errOut << "FATAL error, could not open file pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<"\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    if (genome.chrNameIndex.size()==0) {
        for (uint64 ii=0;ii<genome.nChrReal;ii++) {
            genome.chrNameIndex[genome.chrName[ii]]=ii;
        };
    };

    exonN=0;
    uint64 exonLines=0;
    const vector<vector<string>> exAttrNames({ {genome.pGe.sjdbGTFtagExonParentTranscript}, {genome.pGe.sjdbGTFtagExonParentGene}, genome.pGe.sjdbGTFtagExonParentGeneName, genome.pGe.sjdbGTFtagExonParentGeneType }); //trID, gID, gName, gBiotype
    struct ParsedLine {
        bool exon=false;
        bool chrFound=false;
        bool pastEnd=false;
        string chr1;
        uint64 chr=0, ex1=0, ex2=0;
        char str1=0;
        vector<string> exAttr;
    };

    const size_t blockLines=1<<16;
    vector<string> lines;
    vector<ParsedLine> parsed;
    bool moreLines=true;
    while (moreLines) {
        lines.clear();
        while (lines.size()<blockLines) {
            string oneLine;
            try {
                if (!(genome.personalizedInput ? genome.personalizedInput->nextGtfLine(oneLine)
                                              : static_cast<bool>(getline(sjdbStreamIn,oneLine)))) {
                    moreLines=false;
                    break;
                };
            } catch (const std::exception &error) {
                exitWithError("EXITING because personalized GTF input failed: " + string(error.what()) + "\n",
                              std::cerr,P.inOut->logMain,EXIT_CODE_INPUT_FILES,P);
            };
            lines.push_back(std::move(oneLine));
        };

        parsed.assign(lines.size(),ParsedLine());
        const int parseThreads=genome.personalizedInput ? genome.personalizedInput->availableThreads() : P.runThreadN;
        #pragma omp parallel for schedule(dynamic,1024) num_threads(parseThreads)
        for (long long il=0; il<static_cast<long long>(lines.size()); il++) {
            ParsedLine &line=parsed[static_cast<size_t>(il)];
            string chr1,ddd2,featureType;
            istringstream oneLineStream (lines[static_cast<size_t>(il)]);

            getline(oneLineStream,chr1,'\t');
            getline(oneLineStream,ddd2,'\t');
            getline(oneLineStream,featureType,'\t');
            if (!(chr1.substr(0,1)!="#" && featureType==genome.pGe.sjdbGTFfeatureExon))
                continue;
            line.exon=true;

            if (genome.pGe.sjdbGTFchrPrefix!="-") chr1=genome.pGe.sjdbGTFchrPrefix + chr1;
            line.chr1=chr1;
            const auto contig=genome.chrNameIndex.find(chr1);
            if (contig==genome.chrNameIndex.end())
                continue;
            line.chrFound=true;
            line.chr=contig->second;

            oneLineStream >> line.ex1 >> line.ex2 >> ddd2 >> line.str1 >> ddd2; //read all fields except the last
            if (line.ex2 > genome.chrLength[line.chr]) {
                line.pastEnd=true;
                continue;
            };

            string oneLine1;
            getline(oneLineStream, oneLine1);//get the last field
            replace(oneLine1.begin(),oneLine1.end(),';',' ');//to separate attributes
            replace(oneLine1.begin(),oneLine1.end(),'=',' ');//for GFF3 processing
            replace(oneLine1.begin(),oneLine1.end(),'\t',' ');//replace tabs
            replace(oneLine1.begin(),oneLine1.end(),'\"',' ');//now the only separator is space

            line.exAttr.resize(exAttrNames.size());
            for (uint32 ii=0; ii<exAttrNames.size(); ii++) {
                for (const auto &attr1 : exAttrNames[ii]) {//scan through possible names
                    size_t pos1=oneLine1.find(" " + attr1 + " "); //attribute name is separated by spaces
                    if (pos1!=string::npos)
                        pos1=oneLine1.find_first_not_of(" ", pos1+attr1.size()+1);
                    if (pos1!=string::npos) {
                        line.exAttr[ii]=oneLine1.substr(pos1, oneLine1.find_first_of(" ",pos1)-pos1);
                    };
                };
            };
        };

        for (size_t il=0; il<lines.size(); il++) {
            const string &oneLine=lines[il];
            ParsedLine &line=parsed[il];
            if (!line.exon)
                continue;
            ++exonLines;

            const string &chr1=line.chr1;
            if (!line.chrFound) {//chr not in Genome
                P.inOut->logMain << "WARNING: while processing sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": chromosome '"<<chr1<<"' not found in Genome fasta files for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
                continue; //do not process exons/transcripts on missing chromosomes
            };
            if (line.pastEnd) {
                warningMessage("while processing sjdbGTFfile=" + genome.pGe.sjdbGTFfile + ", line:\n" + oneLine + "\n exon end = " + to_string(line.ex2) +  \
                               " is larger than the chromosome " + chr1 + " length = " + to_string(genome.chrLength[line.chr] ) + " , will skip this exon\n", \
                               std::cerr, P.inOut->logMain, P);
                continue;
            };

            vector<string> &exAttr=line.exAttr;
            if (exAttr[0]=="") {//no transcript ID
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": no transcript_id for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
                exAttr[0]="tr_" + chr1 +"_"+ to_string(line.ex1) +"_"+ to_string(line.ex2) +"_"+ to_string(exonN); //unique name for the transcript
            };

            if (exAttr[1]=="") {//no gene ID
                P.inOut->logMain << "WARNING: while processing pGe.sjdbGTFfile=" << genome.pGe.sjdbGTFfile <<": no gene_id for line:\n";
                P.inOut->logMain << oneLine <<"\n"<<flush;
                exAttr[1]="MissingGeneID";
            };

            if (exAttr[2]=="") {//no gene name
                exAttr[2]=exAttr[1];
            };

            if (exAttr[3]=="") {//no gene name
                exAttr[3]="MissingGeneType";
            };

            transcriptIDnumber.insert(std::pair <string,uint64> (exAttr[0],(uint64) transcriptIDnumber.size()));//insert new element if necessary with a new numeric value
            if (transcriptID.size() < transcriptIDnumber.size()) {//new transcript
                transcriptID.push_back(exAttr[0]);
                if (line.str1=='+') {
                   transcriptStrand.push_back(1);
                } else if (line.str1=='-') {
                   transcriptStrand.push_back(2);
                } else {
                   transcriptStrand.push_back(0);
                };
            };

            geneIDnumber.insert(std::pair <string,uint64> (exAttr[1],(uint64) geneIDnumber.size()));//insert new element if necessary with a $
            if (geneID.size() < geneIDnumber.size()) {//new gene is added
                geneID.push_back(exAttr[1]);
                geneAttr.push_back({exAttr[2],exAttr[3]});
            };

            exonLoci.push_back({ transcriptIDnumber[exAttr[0]],
                                 line.ex1+genome.chrStart[line.chr]-1,
                                 line.ex2+genome.chrStart[line.chr]-1,
                                 geneIDnumber[exAttr[1]] });
            ++exonN;
        };
    };

    if (exonLines==0) {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no ""exon"" lines in the GTF file: " << genome.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file, it must contain some lines with ""exon"" in the 3rd column.\n";
        errOut << "          Make sure the GTF file is unzipped.\n";
        errOut << "          If exons are marked with a different word, use --sjdbGTFfeatureExon .\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };

    if (exonN==0) {
        ostringstream errOut;
        errOut << "Fatal INPUT FILE error, no valid ""exon"" lines in the GTF file: " << genome.pGe.sjdbGTFfile <<"\n";
        errOut << "Solution: check the formatting of the GTF file. One likely cause is the difference in chromosome naming between GTF and FASTA file.\n";
        exitWithError(errOut.str(),std::cerr, P.inOut->logMain, EXIT_CODE_INPUT_FILES, P);
    };
    
    return;
};
