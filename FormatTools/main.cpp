#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include "uRegion.h"
#include "uFunctions.h"
#include "PEComplete.h"
#include "SamTrim.h"
#include "SamSubset.h"


using namespace std;


 /**< Pretty odl..do not use this */


void printHelp();
void completePeSam(string pathname, ostream& output);

int main(int argc, char* argv[])
{

    string firstArg;
    uTagsExperiment expTag;
    ifstream samFile;
    ofstream outputA;
    if (argc <=1){

        cerr << "Arguments needed. -h for help";
        return 0;
    }

    firstArg=argv[1];

    //Case would be wonderful here, but hey, this is C++ and we don't feel like messing with enum and maps.

    if ( firstArg=="-h")
    {
            printHelp();
    }
    else
    if (firstArg=="PEComplete"){
     if(argc <=2)
            {
                cerr<<"Program signature is <filepath>"<<endl;
                return 0;
            }
            string pathname = argv[2];
            completePeSam(pathname, cout);
    }else
    if (firstArg=="TrimSam"){

        if(argc <=4 )
             {
                cerr<<"Program signature is <filepath> <trim left> <trim right>"<<endl;
                return 0;
            }

        string pathname = argv[2];
        int left = atoi(argv[3]);
        int right = atoi(argv[4]);

        SamTrim(left,right, pathname, cout);
    }else
     if (firstArg=="SamToBed"){
     if(argc <=2)
            {
                cerr<<"Program signature is <filepath>"<<endl;
                return 0;
            }
            string pathname = argv[2];
            samFile.open(pathname.c_str());
            cerr <<"Trying to load " << argv[2] << endl;
            expTag.loadFromSam(samFile, true);
            expTag.writeToBed(cout);


          //  getSamSubset();
    }else
     if (firstArg=="SamSubset"){
     if(argc <=2)
            {
                cerr<<"Program signature is <filepath>"<<endl;
                return 0;
            }
            string pathname = argv[2];
            samFile.open(pathname.c_str());
            cerr <<"Trying to load " << argv[2] << endl;
            expTag.loadFromSam(samFile, false);
            string chr;
            int start, end;
            bool stop=false;
            uTagsChrom tempChrom;
            while(!stop){
              cerr << "Region chr start and end. type stop to end program" << endl;
              cin >> chr;

              if (chr=="stop"){
                stop=true;
                break;
              }
              cin >> start;
              cin >> end;

            tempChrom= expTag.getSubset(chr,start,end);
            tempChrom.writeSamToOutput(cout);
            }

          //  getSamSubset();
    }else
    if (firstArg=="GetSignal"){

     if (argc <=3){
        cerr << "Arguments are <inputSam> <outputFile> ";
        return 0;
    }
        samFile.open(argv[2]);
        outputA.open(argv[3]);

        cerr <<"Trying to load " << argv[2] << endl;
        expTag.loadFromSam(samFile, true);

        string chr;
        int start, end;
        bool bComment=false;;
        vector<float> Signal;
        uRegion ourRegion;
        cout << "Type Y if you want a comment line for each signal" << endl;
        cin>>chr;
        if (chr.at(0)=='Y')
            bComment=true;
        while(true){
              cout << "Region chr start and end" << endl;
              cin >> chr;
              cin >> start;
              cin >> end;

              ourRegion.setChr(chr);
              ourRegion.setStart(start);
              ourRegion.setEnd(end);

              Signal.clear();
              Signal.resize(ourRegion.getLenght());
              RegionTags::getSignal(&ourRegion, &expTag, true);
              Signal=ourRegion.getSignal();
              if (bComment)
                    outputA<< "Signal from " << chr << " " << start << " " << end << endl;

              for (int i=0; i< (int)Signal.size(); i++)
                 outputA << Signal.at(i) << "\t";

                outputA<<endl;
         }
   }

}

void printHelp(){
    cerr <<"PEComplete"<<endl;
    cerr <<"TrimSam "<<endl;
    cerr <<"GetSignal"<<endl;
    cerr <<"SamToBed"<<endl;
    cerr <<"SamSubset"<<endl;
}
