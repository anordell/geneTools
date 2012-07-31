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

int findNearestTag( uChipPeak ourTag,int pos, VecChipPeak* ourTagsinPeak);

int main(int argc, char* argv[])
{
    string bedFileName;
    string peakFileName;
    string outputFileName;
    ifstream inputBed,inputPeaks;
    ofstream outputBed;
    vector<PeakTags> vecPeaktags;
    PeakTags tempPeak;

  //  if (argc < 4){
  //      cout<< " Signature is  <bedFile> <Peakfile> <OutputName> ";
       // return 0;
  //  bedFileName=argv[1];
  //  peakFileName=argv[2];
    outputFileName="resulstsall.txt";

    bedFileName ="H2AZq10cutsam.bed";
  //  peakFileName="H2AZq10_withdup.bed";


    ChipPeakVectorMap mymap;

    uGeneStats geneStats,geneStats2;
    uChipSeq tagdata, peakdata;
    uGeneRef geneRef;

    inputBed.open(bedFileName.c_str());
    outputBed.open(outputFileName.c_str());
  //  inputPeaks.open(peakFileName.c_str());

  //  outputBed.open(outputFileName.c_str());

     cout << "Start loading \n";

    tagdata.loadFromBed(inputBed);

    cout << "Finished loading \n";

//    tagdata.outputBedFormat(cout);


    VecChipPeak* pVecTags;
    ChipPeakVectorMap ourMap, tagMap;
    ChipPeakVectorMap::iterator iterMap;

    vector <int >distance;
    int totaldist,curdist;

    ourMap= tagdata.getMap();
    //For every chromosome
    cout << "started Chrom \n";
     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
         cout << "started Chrom"<< iterMap->first<< " \n";
        pVecTags=(&iterMap->second);
       //For every tag on that chromosome
        for (unsigned int i=0; i< pVecTags->size();i++){
            curdist=findNearestTag( pVecTags->at(i),i,  pVecTags);
            distance.push_back(curdist);
            totaldist+=curdist;
         }

    }
    totaldist=totaldist/ tagdata.count();


       cout << "Finished Chrom \n";
    int groupedistant[10000];

    for (int i=0; i< 10000; i++){

       groupedistant[i]=0;


     }

    for (int i=0; i< distance.size(); i++){

        if (distance.at(i)>=10000)
            groupedistant[10000]++;
         else
            groupedistant[distance.at(i)]++;
        }


cout << "Writing \n";
     for (int i=0; i<10000; i++){

        outputBed << "Group " << i << "\t" ;
        outputBed << groupedistant[i] <<"\n";


     }


   /* for (int i=0; i< distance.size(); i++){

        cout << "Tag " << i << "\t" ;
        cout << distance.at(i) <<"\n";
    } */

    return 0;
}

int findNearestTag( uChipPeak ourTag,int pos, VecChipPeak* ourTagsinPeak){

    int distance = 999999;
    int distprev=0, disnext=0;
    int prev,next;
    bool foundself=false;
    prev= pos;
    next=pos;

    //find previous
    bool pFound=false, nFound=false;
    while (pFound==false){
        prev--;
        if (prev!=-1)
        {
             if(ourTag.getchromStart()!= ourTagsinPeak->at(prev).getchromStart() ){
                distprev= abs(ourTag.getchromStart()-ourTagsinPeak->at(prev).getchromStart());
                 pFound=true;
             }
        }
        else{
            pFound=true;
            distprev=999999;
        }
    }

  while (nFound==false){
      next++;
        if (next< ourTagsinPeak->size())
        {
             if(ourTag.getchromStart()!= ourTagsinPeak->at(next).getchromStart() ){
                disnext= abs(ourTag.getchromStart()-ourTagsinPeak->at(next).getchromStart());
                 nFound=true;
             }
        }
        else{
            nFound=true;
            disnext=999999;
        }
    }


    if (disnext < distprev )
        distance = disnext;
    else
        distance = distprev;


/*
     for (int i=0; i< ourTagsinPeak.size(); i++){

        tempdist= abs(ourTag.getchromStart()-ourTagsinPeak.at(i).getchromStart());

        if ( ourTag.getchromStart()== ourTagsinPeak.at(i).getchromStart())
           foundself=true;

        if (tempdist< distance){
            if( ourTag.getchromStart()!= ourTagsinPeak.at(i).getchromStart() )
                distance = tempdist;
        }

            if (foundself==true){
                if (tempdist!=0)
                    break;
            }
    }
*/
    return distance;
}
