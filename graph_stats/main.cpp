#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>
#include "uGeneStats.h"

using namespace std;

int main(int argc, char* argv[])
{

    ifstream inputBed,inputBed2,inputGenes;
    string bedFileName;

    ChipPeakVectorMap mymap;
    vector<int> windowResult;
      bedFileName=argv[1];

    uGeneStats geneStats,geneStats2;
    uChipSeq chipTest,chipTest2;
    uGeneRef geneRef;

    inputBed.open(bedFileName.c_str());
 //   inputBed2.open("C:\\Work\\data\\ours\\H2AZ\\H2AZ\\peakcalling\\H2AZ-2_peaks.bed");
    inputGenes.open("C:\\Work\\tools\\junkyard\\junktesting\\testdata\\UCSCCanonAll.bed");


    chipTest.loadFromBed(inputBed);
   // chipTest2.loadFromBed(inputBed2);
    geneRef.loadfromUCSCDB(inputGenes);
    geneStats.loadData(geneRef,chipTest);
  //  geneStats2.loadData(geneRef,chipTest2);


/*
    geneRef.loadfromUCSCDB(inputBed);
    chipTest.loadFromBed(inputBed2);

    geneStats.loadData(geneRef, chipTest);

    int count, countGene;
    count= geneStats.countTss(300);
    countGene= geneStats.countInGene();

    cout <<"Number of peaks:"<<chipTest.count()<<"\n";
    cout <<"Number of peaks at tss" <<count<<"\n";
    cout <<"Number of peaks in Gene" <<countGene<<"\n";

*/
/*
    cout << "H2AZE2\n";
    cout <<"Peaks Count "<< chipTest.count()<< "\n";
    cout <<"Peaks overlapping TSS "<< geneStats.countOverlapTss()<<"\n";
    cout <<"Middle peaks in 1KP of TSS"<< geneStats.countTss(1000)<<"\n";
    cout <<"Middle peaks in gene "<< geneStats.countInGene()<<"\n";
    cout <<"Middle peaks in exon" <<geneStats.countInExon()<<"\n";
    cout <<"Middle peaks in Proximal promoter of 1K size "<<geneStats.countinProximalProm(1000)<<"\n";
    cout <<"Middle peaks in Proximal promoter of 3K size "<<geneStats.countinProximalProm(3000)<<"\n";


    cout << "H2AZ\n";
    cout <<"Peaks Count "<< chipTest2.count()<< "\n";
    cout <<"Peaks overlapping TSS "<< geneStats2.countOverlapTss()<<"\n";
     cout <<"Middle peaks in 1KP of TSS"<< geneStats2.countTss(1000)<<"\n";
    cout <<"Middle peaks in gene "<< geneStats2.countInGene()<<"\n";
    cout <<"Middle peaks in exon" <<geneStats2.countInExon()<<"\n";
    cout <<"Middle peaks in Proximal promoter of 1K size "<<geneStats2.countinProximalProm(1000)<<"\n";
     cout <<"Middle peaks in Proximal promoter of 3K size "<<geneStats2.countinProximalProm(3000)<<"\n";
    // cout <<"Unique tags on CTCF-CD4:"<<chipTest.countUnique()<<"\n";
    // cout <<"Unique tags on H2AZme1-E2:"<<chipTest2.countUnique()<<"\n";
*/

    windowResult=geneStats.getTssProximityWindow();

    cout <<bedFileName<< "\n";
    for(unsigned int k=0; k< windowResult.size(); k++)
        cout << windowResult.at(k)<<"\t";

     cout << "\n";
    return 0;
}
