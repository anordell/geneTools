#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>
#include "uGeneStats.h"

using namespace std;




struct PeakTags{
    uChipPeak ourPeak;
    vector< uChipPeak> ourTagsinPeak;
};



int findNearestTag( uChipPeak ourTag, vector< uChipPeak> ourTagsinPeak, int myPos);
int avgTagDensity( PeakTags ourPeakTag);

int main(int argc, char* argv[])
{
    string bedFileName;
    string peakFileName;
    string outputFileName;
    ifstream inputBed,inputPeaks;
    ofstream outputBed;
    vector<PeakTags> vecPeaktags;
    PeakTags tempPeak;

    if (argc < 4){
  //      cout<< " Signature is  <bedFile> <Peakfile> <OutputName> ";
       // return 0;
    }
  //  bedFileName=argv[1];
  //  peakFileName=argv[2];
  //  outputFileName=argv[3];
    int totalcount=0;
    int curcoun=0;

    bedFileName ="H2AZtagsoverlap.bed";
    peakFileName="H2AZq10_withdup.bed";


    ChipPeakVectorMap mymap;

    uGeneStats geneStats,geneStats2;
    uChipSeq tagdata, peakdata;
    uGeneRef geneRef;

    inputBed.open(bedFileName.c_str());

    inputPeaks.open(peakFileName.c_str());

    outputBed.open(outputFileName.c_str());

    tagdata.loadFromBed(inputBed);
    peakdata.loadFromBed(inputPeaks);


    VecChipPeak* pVecPeak;
    VecChipPeak* pVecTags;
    ChipPeakVectorMap ourMap, tagMap;
    ChipPeakVectorMap::iterator iterMap;

    tagMap= tagdata.getMap();
    ourMap= peakdata.getMap();
    //For every chromosome
     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
        pVecPeak=(&iterMap->second);
       //For every peak on that chromosome
        for (unsigned int i=0; i< pVecPeak->size();i++){
            tempPeak.ourPeak=pVecPeak->at(i);
            tempPeak.ourTagsinPeak.clear();

            //Get chromosome

            //For every tag on our chromosome, is it within our peak, if so include
            pVecTags= &tagMap[iterMap->first];

                for (int x=0; x<pVecTags->size(); x++ ){

                    if ( checkOverlap(tempPeak.ourPeak.getchromStart(),tempPeak.ourPeak.getchromEnd(),pVecTags->at(x).getchromStart(),pVecTags->at(x).getchromEnd() ) )
                    {
                        int start, end;
                        start= pVecTags->at(x).getchromStart();
                        end = pVecTags->at(x).getchromEnd();
                        tempPeak.ourTagsinPeak.push_back(pVecTags->at(x));
                    }
                }

         vecPeaktags.push_back(tempPeak);

         }

    }

    for (int i=0; i< vecPeaktags.size(); i++){

        cout << "Peak " << i << "\t" ;
        cout << avgTagDensity(vecPeaktags.at(i)) <<"\n";
    }

    return 0;
}

//Measure tag density of our peak.
int avgTagDensity( PeakTags ourPeakTag){

    int density=0;

  vector< uChipPeak> tempVec;
    for (int i=0; i< ourPeakTag.ourTagsinPeak.size(); i++){

       // tempVec=ourPeakTag.ourTagsinPeak;
       // tempVec.erase(tempVec.begin()+i);

    density+=  findNearestTag(ourPeakTag.ourTagsinPeak.at(i), ourPeakTag.ourTagsinPeak, i);



    }
    density = density/ourPeakTag.ourTagsinPeak.size();

return density;
}

int findNearestTag( uChipPeak ourTag, vector< uChipPeak> ourTagsinPeak, int myPos){

    int distance = 20000;
    int tempdist;

     for (int i=0; i< ourTagsinPeak.size(); i++){

         tempdist= abs(ourTag.getchromStart()-ourTagsinPeak.at(i).getchromStart());
        if (tempdist< distance){
            if( (i!= myPos)&&( ourTag.getchromStart()!= ourTagsinPeak.at(i).getchromStart() ) )
                distance = tempdist;

        }

    }

    return distance;
}




