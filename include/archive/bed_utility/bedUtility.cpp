
#include "bedFile.h"
#include "bedUtility.h"
#include "chromConstants.h"

using namespace std;

bedFile overlapBedfiles(bedFile* fileA, bedFile* fileB, int overlapcount,int distance, int overlap_para){

//Return the elements of fileA that overlap an element in fileB.
//Partial (0) returns only the part of th element overlapping
//Complete (1) returns the entire element.

//The bedFile we will return

bedFile overlapResults;
TagRead currentTagA, currentTagB,tempRead;
int start,end,size;
bool found=false;

//For every chromosome
 for (int i=0; i<CHROMCOUNT; i++)
     {
            ChromSum *pChromA, *pChromB;
            pChromA = fileA->getChrom(i);
            pChromB = fileB->getChrom(i);
            //Every if elements form file A overlap with
           for (unsigned int k=0; k <pChromA->tList.size();k++)
           {
                found=false;
                currentTagA = pChromA->tList.at(k);
                //Spread tag if for distance

                //Elements from fileB
                for (unsigned int n=0; n <pChromB->tList.size();n++)
                {
                    currentTagB =   pChromB->tList.at(n);

                    if ( (( currentTagA.chromStart-distance) < ( currentTagB.chromEnd) ) &&
                      ((currentTagA.chromEnd+distance) > ( currentTagB.chromStart)    )  )
                      {
                          //We start counting at the farthest and stop at the nearest
                          if ((currentTagA.chromStart-distance)>currentTagB.chromStart)
                                start= (currentTagA.chromStart-distance);
                          else
                                start =currentTagB.chromStart;
                          if ((currentTagA.chromEnd+distance)<currentTagB.chromEnd)
                                end= (currentTagA.chromEnd+distance);
                          else
                                end =currentTagB.chromEnd;

                        size = (end-start);
                            //If not overlap count, we don't add
                        if (size>=overlapcount)
                        {
                            found=true;
                            //Return only the overllaping portion


                            switch (overlap_para){
                                //return the overlapping portion including flanks
                                case RETURN_PARTIAL:
                                    tempRead.chrom=i;
                                    tempRead.chromStart=start;
                                    tempRead.chromEnd=end;
                                    break;
                            //Return only the original site
                                case RETURN_COMPLETE:
                                    tempRead= currentTagA;
                                    break;

                                    //return the site plus flanks.
                                case RETURN_COMPLETE_FLANK:
                                    tempRead.chrom=i;
                                    tempRead.chromStart=(start-distance);
                                    if (tempRead.chromStart<0)
                                        tempRead.chromStart=0;

                                    tempRead.chromEnd=(end+distance);
                                    if (tempRead.chromEnd>CHROMSIZE[i])
                                        tempRead.chromEnd=CHROMSIZE[i];

                                    break;
                            }

                            overlapResults.addReadtoChrom(tempRead, i);
                        }

                      }
                  //if we found
                if (found== true)
                    break;
                }
            //End TagsA of chromosome
            }

   //End Chromosome
    }



 return overlapResults;

}



bedFile substractBedfiles(bedFile* fileA, bedFile* fileB, int overlapcount,int distance){

//Return the elements of fileA that do not overlap an element in fileB
//The bedFile we will return
bedFile substractResults;
TagRead currentTagA, currentTagB,tempRead;
bool found=false;

//For every chromosome
 for (int i=0; i<CHROMCOUNT; i++)
     {
            ChromSum *pChromA, *pChromB;
            pChromA = fileA->getChrom(i);
            pChromB = fileB->getChrom(i);
            //Every if elements form file A overlap with
           for (unsigned int k=0; k <pChromA->tList.size();k++)
           {
                found=false;
                currentTagA = pChromA->tList.at(k);
                //Spread tag if for distance

                //Elements from fileB
                for (unsigned int n=0; n <pChromB->tList.size();n++)
                {
                   currentTagB =   pChromB->tList.at(n);
                   tempRead=isoverlapTag( currentTagA, currentTagB, overlapcount, distance, RETURN_COMPLETE);

                   if (tempRead.chromEnd){
                       found=true;
                       break;
                   }
                }
            //We return TagA
                if (found==false)
                    substractResults.addReadtoChrom(currentTagA, i);

                found=false;
            //End TagsA of chromosome
            }
   //End Chromosome
    }
 return substractResults;
}


//Return values
//0 Not
//1 True
//More to come maybe
TagRead isoverlapTag(TagRead tagA, TagRead tagB, int overlapcount, int distance,int overlap_para)
{
    int start,end,size;
    bool found=false;
    TagRead tempRead;

     if ( (( tagA.chromStart-distance) <= ( tagB.chromEnd) ) &&
      ((tagA.chromEnd+distance) >= ( tagB.chromStart)    )  )
      {
          //We start counting at the farthest and stop at the nearest
          if ((tagA.chromStart-distance)>tagB.chromStart)
                start= (tagA.chromStart-distance);
          else
                start =tagB.chromStart;
          if ((tagA.chromEnd+distance)<tagB.chromEnd)
                end= (tagA.chromEnd+distance);
          else
                end =tagB.chromEnd;

        size = (end-start);
            //If not overlap count, we don't add
        if (size>=overlapcount)
        {
            found=true;
            //Return only the overllaping portion
            if ( overlap_para== RETURN_PARTIAL)
            {
                tempRead.chrom=tagA.chrom;
                tempRead.chromStart=start;
                tempRead.chromEnd=end;
            }
            else //Return all
            {
                tempRead= tagA;
            }
        }
      }
    if (found==false){
        tempRead.chrom=0;
        tempRead.chromStart=0;
        tempRead.chromEnd=0;
    }

return tempRead;

}


barrierResults barrierOverlap(bedFile inputData, bedFile comparisonData, int distance){

barrierResults resultData;



//For every chromosomes
for(int i=0; i< CHROMCOUNT;i++ ){

    //Verify every tag in the comparisonData to everyislands in the inputData

    for (unsigned int k=0; k<comparisonData.getChrom(i)->tList.size();k++){
        bool found;
        TagRead currentTag;
        found=false;

        currentTag = comparisonData.getChrom(i)->tList[k];


        for(unsigned int n=0; n<inputData.getChrom(i)->iList.size(); n++ ){
            TagIsland currentIslands;
            currentIslands= inputData.getChrom(i)->iList[n];
            //If tags overlap, we check, is it in islands or a boundary?
            if ( (currentTag.chromStart < ( currentIslands.iStop +distance) ) &&
              (currentTag.chromEnd > ( currentIslands.iStart -distance)    )  )
              {
                  //check if central
                  found = true;
                  if ((currentTag.chromStart > (currentIslands.iStart+distance) )&&
                     (currentTag.chromEnd < (currentIslands.iStop-distance) ) ){

                        resultData.centerTags[i].push_back(currentTag);
                     }
                     //POtherwise barriers
                     else{
                          resultData.barrierTags[i].push_back(currentTag);
                     }
              }
              //If not normal tag, exist for.
        if (found== true)
            break;

        }

        //Otherwise put in list of non-barriers tags.
        if (found==false)
            resultData.otherTags[i].push_back(currentTag);
    }
}

return resultData;

}



string getChromfromint(int chromNum){


    string outputChrom;

    if (chromNum>=CHROMX_POS){
               if (chromNum== CHROMX_POS)
                    outputChrom="chrX";
                if (chromNum== CHROMY_POS )
                    outputChrom="chrY";
                if (chromNum== CHROMM_POS )
                    outputChrom="chrM";

           }
           else{
               outputChrom= "chr"+convertInt(chromNum+1);

           }
    return outputChrom;
}


int getChromfromString(string chromName){

    int chromNum;
    //Ignore the first three letters
   // cout << chromName;
    chromName = chromName.substr(3);
    chromNum  = atoi(chromName.c_str());
    if ((chromNum==0)||(chromName.length()>4)){
        if (chromName=="M")
            chromNum=CHROMM_POS;
            else
            if (chromName=="X")
                chromNum=CHROMX_POS;
                else
                if (chromName=="Y")
                    chromNum=CHROMY_POS;
                    else
                    return -1;

     }
     else{
    chromNum--; }

    return chromNum;
}
