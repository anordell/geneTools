#include <iostream>
#include <istream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>
#include "uChipSeq.h"

using namespace std;


float percentBig(int A, int B);
//We do Biggest over smallest And output size differential.
int main(int argc, char* argv[])
{
    uChipSeq chipSeqA, chipSeqB;
    string bedFileNameA, bedFileNameB;
    ifstream inputBedA,inputBedB;
    int extend=0;
     if (argc < 3){
        cout<< " Signature is  <bedFileA> <BedFileB>";
        return 0;
    }

    bedFileNameA=argv[1];
    bedFileNameB=argv[2];

    if (argc==4)
        extend= atoi(argv[3]);


   inputBedA.open(bedFileNameA.c_str());
   inputBedB.open(bedFileNameB.c_str());




//If fail opening.
if  ((inputBedA.fail())||(inputBedB.fail()))  {


          cout << "Error opening files.\n";
          cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    }else{
        float percentOverlapAB, percentOverlapBA;
         chipSeqA.loadFromBed( inputBedA);
        chipSeqB.loadFromBed( inputBedB);

        //Smallest on Biggest

            percentOverlapBA=chipSeqB.CorPeaks(chipSeqA,extend);

            percentOverlapAB= chipSeqA.CorPeaks(chipSeqB,extend);



        cout << percentOverlapAB<<"\n";
        cout << percentOverlapBA<<"\n";
        cout << chipSeqA.count()<<"\n";
        cout << chipSeqB.count()<<"\n";
        cout << percentBig( chipSeqA.count(),chipSeqB.count())<<"\n";
    }


    return 0;
}


float percentBig(int A, int B){

    float result;

    if (A>B)
        result = (float)A/(float) B;
    else
        result = (float)B/(float)A;


return result;
}
