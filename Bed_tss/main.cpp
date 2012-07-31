#include <bedFile.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

  // testBed.mapTagstoWindows(1000);
  //  testBed.outputWindows("CD4-H3K27me3_1000_windows_beta.bed");

   // testBed.countIslands(35,4,7);
    //testBed.outputIslands("CD4-H3K27me3_35_4_7_islands_beta2.bed");


int main(int argc, char* argv[])
{
 bedFile bedA,bedB, resultBed;
int parameter;

bedFile tssTransform(bedFile bedData);
bedFile StartEndTransform(bedFile);

string fileA, fileOutput;

fileA=argv[1];
fileOutput=argv[2];
parameter=atoi(argv[3]);

//fileA="h3k27barcomp.txt";
//fileOutput="happy.txt";
//parameter=1;

argc=4;

if (argc<3){

        cout << "Usage: Bed_tss <filename> <1:0> \n";
        cout << "0 for start. 1 for start and end transform \n";

        return 0;
}

    if (bedA.readTags(fileA))  {

        bedA.mapTagstoChrom();
        bedA.sortChromTags();

    if (parameter=0)
        bedA= tssTransform(bedA);
    if (parameter=1)
         bedA= StartEndTransform(bedA);
    }


    bedA.outputBedTags(fileOutput);

    return 0;
}


bedFile tssTransform(bedFile bedData)

{

       for (int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        pCurChrom = bedData.getChrom(i);
         for (unsigned int k=0; k <pCurChrom->tList.size();k++){
                    pCurChrom->tList.at(k).chromEnd=(pCurChrom->tList.at(k).chromStart+4 );
            }

        }
    return bedData;
}

bedFile StartEndTransform(bedFile bedData)

{
TagRead tempread;
bedFile tempBed;

    for (int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        pCurChrom = bedData.getChrom(i);
         for (unsigned int k=0; k <pCurChrom->tList.size();k++){

             //Add start
             tempread= pCurChrom->tList.at(k);
             tempread.chromEnd=tempread.chromStart+3;

             tempBed.addReadtoChrom(tempread,i);

             //End
             tempread=pCurChrom->tList.at(k);
             tempread.chromStart=tempread.chromEnd-3;

            tempBed.addReadtoChrom(tempread,i);

            }
        }

return tempBed;

}

