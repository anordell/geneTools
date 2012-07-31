#include "uGeneStats.h"

uGeneStats::uGeneStats(uGeneRef ourRef, uChipSeq ourChip) {  // default constructor

 ourGeneRef= ourRef;

 ourChipSeq= ourChip;

}

uGeneStats::uGeneStats() {  // default constructor


}

uGeneStats::~uGeneStats() {  // default destructor

}

void uGeneStats::loadData(uGeneRef ourRef, uChipSeq ourChip){


 ourGeneRef= ourRef;

 ourChipSeq= ourChip;


}

int uGeneStats::countTss(int size){


    VecChipPeak* pVecPeak;
    int count =0;
    int pos;
    ChipPeakVectorMap ourMap;
    ChipPeakVectorMap::iterator iterMap;

    ourMap =ourChipSeq.getMap();

     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
        pVecPeak=(&iterMap->second);
        //cout << &iterMap->first << " : " << pVecPeak->size() << "\n";
        for (unsigned int i=0; i< pVecPeak->size();i++){
           pos = ((pVecPeak->at(i).getchromStart()+pVecPeak->at(i).getchromEnd())/2);
           count+= ourGeneRef.isTss(iterMap->first,pos,size);
        }
    }

return count;
}

 int uGeneStats::countInGene(){

   VecChipPeak* pVecPeak;
    int count =0;
    int pos;
    ChipPeakVectorMap ourMap;
    ChipPeakVectorMap::iterator iterMap;

    ourMap =ourChipSeq.getMap();

     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
        pVecPeak=(&iterMap->second);
        //cout << &iterMap->first << " : " << pVecPeak->size() << "\n";
        for (unsigned int i=0; i< pVecPeak->size();i++){
           pos = ((pVecPeak->at(i).getchromStart()+pVecPeak->at(i).getchromEnd())/2);
           count+= ourGeneRef.isGene(iterMap->first,pos);
        }
    }


return count;
}

int uGeneStats::countinProximalProm(int size){

 VecChipPeak* pVecPeak;
    int count =0;
    int pos;
    ChipPeakVectorMap ourMap;
    ChipPeakVectorMap::iterator iterMap;

    ourMap =ourChipSeq.getMap();

     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
        pVecPeak=(&iterMap->second);
        //cout << &iterMap->first << " : " << pVecPeak->size() << "\n";
        for (unsigned int i=0; i< pVecPeak->size();i++){
           pos = ((pVecPeak->at(i).getchromStart()+pVecPeak->at(i).getchromEnd())/2);
           count+= ourGeneRef.isProximalProm (iterMap->first,pos,size);
        }
    }

return count;

}

int uGeneStats:: countInExon(){

VecChipPeak* pVecPeak;
    int count =0;
    int pos;
    ChipPeakVectorMap ourMap;
    ChipPeakVectorMap::iterator iterMap;

    ourMap =ourChipSeq.getMap();

     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
        pVecPeak=(&iterMap->second);
        //cout << &iterMap->first << " : " << pVecPeak->size() << "\n";
        for (unsigned int i=0; i< pVecPeak->size();i++){
           pos = ((pVecPeak->at(i).getchromStart()+pVecPeak->at(i).getchromEnd())/2);
           count+= ourGeneRef.isExon(iterMap->first,pos);
        }
    }

return count;
}
//Count how many sites overlaptheTSS
int uGeneStats::countOverlapTss(){

 VecChipPeak* pVecPeak;
    int count =0;
    ChipPeakVectorMap ourMap;
    ChipPeakVectorMap::iterator iterMap;

    ourMap =ourChipSeq.getMap();

     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
        pVecPeak=(&iterMap->second);
        for (unsigned int i=0; i< pVecPeak->size();i++){
           //pos = ((pVecPeak->at(i).getchromStart()+pVecPeak->at(i).getchromEnd())/2);
           count+= ourGeneRef.overlapsTss(iterMap->first,pVecPeak->at(i).getchromStart(),pVecPeak->at(i).getchromEnd());
        }
    }

return count;



}

//Doto
 //Count how many peaks overlap a Gene
int uGeneStats::countOverlapgene(){





return 0;
    }

//Doto
vector<uGene> getNearestGenes(int distance){

vector <uGene> resultvector;


return resultvector;
}

//Return a vector containing windows of 5bp each and the number of peaks that overlap that window.

vector<int> uGeneStats::getTssProximityWindow(){

    VecChipPeak* pVecPeak;

    //Vector is 3K before tss and 1 K after. So 4000/5= 800 buckets, with 0 from -3000 to 2955 and 799 from 995 to 1000
    vector<int> proxWindow;

    proxWindow.resize(800);

    int distTss, startTrans,EndTrans,curTrans;
    ChipPeakVectorMap ourMap;
    ChipPeakVectorMap::iterator iterMap;
    uChipPeak tempPeak;

    ourMap =ourChipSeq.getMap();

    //On every chrom
     for (iterMap = ourMap.begin(); iterMap != ourMap.end(); ++iterMap) {
        pVecPeak=(&iterMap->second);
        //For every peak on the chrom, fill in our bucket.
        for (unsigned int i=0; i< pVecPeak->size();i++){
          tempPeak=pVecPeak->at(i);
          distTss= ourGeneRef.distanceFromClosestTss(iterMap->first, tempPeak.getchromStart());

          //if closest TSS is not within our bucket, we ignore, else we process.
            if ((distTss>1000)||( (distTss+tempPeak.getSize()) < -3000 )){
             //Do nothing
            }
            //Transform our data and add to bucket.
            else{

            startTrans=distTss+3000;
            EndTrans=startTrans+tempPeak.getSize();
            curTrans=startTrans;
            //Distance of 0 is at left promoter edge, etc
            //We walk forward.
            while(curTrans<=EndTrans){
                    if ((curTrans>=0)&&(curTrans<=3999)){

                        proxWindow.at(curTrans/5)++;

                    }
                    curTrans+=5;
                }
            }
        }
    }
return proxWindow;

}
