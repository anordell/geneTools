#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>
#include "uGeneStats.h"

using namespace std;

int main(int argc, char* argv[])
{
    string bedFileName;
    string geneFileName;
    string outputFileName;
    ifstream inputBed,inputGenes;
    ofstream outputBed;


    if (argc < 4){
        cout<< " Signature is  <bedFile> <GenomeSummaryFile> <OutputName> ";
        return 0;
    }
    bedFileName=argv[1];
    geneFileName=argv[2];
    outputFileName=argv[3];


    ChipPeakVectorMap mymap;

    uGeneStats geneStats,geneStats2;
    uChipSeq chipData;
    uGeneRef geneRef;

    inputBed.open(bedFileName.c_str());

    inputGenes.open(geneFileName.c_str());
    outputBed.open(outputFileName.c_str());

    chipData.loadFromBed(inputBed);
    geneRef.loadfromUCSCDB(inputGenes);
    geneStats.loadData(geneRef,chipData);


  /*  outputBed <<"Peaks Count "<< chipData.count()<< "\n";
    outputBed <<"Average peak " << chipData.avgPeakSize()<< "\n";
    outputBed <<"Largest peak " << chipData.maxPeakSize()<< "\n";
    outputBed <<"Smallest peak " << chipData.minPeakSize()<< "\n";
    */
    chipData.printStats(outputBed);
    outputBed << bedFileName<<"\n";;
    outputBed <<"Peaks overlapping TSS"<< "\t"<< geneStats.countOverlapTss()<<"\n";
    outputBed <<"Middle peaks in 1KP of TSS "<< "\t"<<geneStats.countTss(1000)<<"\n";
    outputBed <<"Middle peaks in gene "<< "\t"<<geneStats.countInGene()<<"\n";
    outputBed <<"Middle peaks in exon" <<"\t"<<geneStats.countInExon()<<"\n";
    outputBed <<"Middle peaks in Proximal promoter of 1K size "<< "\t"<<geneStats.countinProximalProm(1000)<<"\n";
    outputBed <<"Middle peaks in Proximal promoter of 3K size "<< "\t"<< geneStats.countinProximalProm(3000)<<"\n";




    return 0;
}
