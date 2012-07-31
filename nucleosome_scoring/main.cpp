#include "uTags.h"
#include "uNucleoBin.h"
#include "clustering.h"
#include "uRegion.h"
#include "scorenorm.h"
#include "uFunctions.h"
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "functions.h"
#include "score.h"
#include "filter.h"
#include "NScore.h"
using namespace std;

//Two main objectives
//Generate the score implement in the Alz Nature paper of 2011
//Decompose our signal and see what it does.
void printHelp();
void test(int argc, char* argv[]);
void outputLocalMaximalBedGraph(string pathname);
void densityCount(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    string firstArg;

   // cerr <<"Roughly "  << ((float)5/(float)4)*100.f << "%"<< endl;

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
    if ( firstArg=="scoreandnorm"){
        scoreAndNormalize(argc, argv);
    }
    else
    if ( firstArg=="score"){
        score(argc, argv);
    }
    else
    if (firstArg=="densityCount"){
        densityCount(argc,argv);
    }
    else
    if ( firstArg=="filter"){
        filter(argc,argv);
    }
    else
    if ( firstArg=="Nscore"){
        NScore(argc,argv);
    }
    else
    if( firstArg=="maxima"){

        if(argc <=2)
        {
            cerr<<"Program signature is <bedGraph> " <<endl;
            return 0;
        }
        string pathname = argv[2];
        outputLocalMaximalBedGraph(pathname);
    }
    else
    if (firstArg=="decomp"){
         decomp(argc, argv);
    }
    else
    if (firstArg=="test"){
         test(argc, argv);
    }
    else{
       cerr <<"Invalid arguments. -h for help" <<endl;;
    }

 return 0;
}
void printHelp(){

    cerr << "Possible argument footprints are:" <<endl;
    cerr << "score -f <filepath> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsComplete=0]  " << endl;
    cerr << "Nscore -f <filepath> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsComplete=0] " << endl;
    cerr << "filter -b <bedGraph> -s [score threshold=0] -t [top threshold=2]" << endl;
    cerr << "decomp -s <samFile> -c <chromosone> -b <BedGraph Maxima> -o <SamLeft> <SamDecomp> ";
    cerr << "maxima <bedGraph>" << endl;
   // cerr << "rec-score <samfile> <chromosomename> <location bedGraph>" << endl;
    cerr << "output is to standard output : Do not forget to  pipe! ( > results.txt )" <<endl;

}
void test(int argc, char* argv[]){

    uNucleoBin test,testloadb;


    ifstream testa("/media/Data/Work/data/nucleosome/E2/specMaps/nuc-mapE2.chr18.Nscore.bedgraph");
    ifstream testb("/media/Data/Work/data/nucleosome/NoE2/specMaps/nuc-mapNoE2.chr18.Nscore.bedgraph");


    test.loadScoreFromBedGraph(testa);
    testloadb.loadScoreFromBedGraph(testb);

    if (test.DensityMapScore.size()>testloadb.DensityMapScore.size())
        test.DensityMapScore.resize(testloadb.DensityMapScore.size());
    if (testloadb.DensityMapScore.size()>test.DensityMapScore.size())
        testloadb.DensityMapScore.resize(test.DensityMapScore.size());
    cout << test.DensityMapScore.size()<<" " << testloadb.DensityMapScore.size() <<endl;
    utility::pause_input();
    auto result= clustering::hauftsmanTwoBins(test,testloadb);

    cout <<test.DensityMapScore.size() <<endl;
    utility::pause_input();

}

//Write the number of counts for each position.
void densityCount(int argc, char* argv[]){

    if(argc <=4)
        {
            cerr<<"Program signature is -s <samfile> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsComplete=0]";
            return;
        }
        int threshold=0;
        bool isComplete=false;
        string pathname, chromName;// = argv[2];
        string bedpath;
    for (int i = 1; i < argc; i++) {
                if (strcmp(argv[i], "-s")==0) {
                    // We know the next argument *should* be the filename:
                    pathname = argv[i + 1];
                } else if (strcmp(argv[i],"-c")==0) {
                    chromName = argv[i + 1];
                }  else if (strcmp(argv[i],"-n")==0) {
                  threshold = atoi(argv[i + 1]);
                }
                else if(strcmp(argv[i],"-p")==0){
                  int temp = atoi(argv[i+1]);
                    if (temp)
                        isComplete=true;
                }
        }
 if ((pathname.size()==0)||(chromName.size()==0)){
        cerr<<"Program signature is -f <filepath> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsPe=0]";
        return;
    }

    uNucleoBin ourBin;
    ourBin.setChromName(chromName);
    cerr << "Trying to Score out file" << endl;
    mainScoring(pathname, ourBin, threshold,isComplete);
    ourBin.writeCounttoBedGraph(cout, 30);

}


//Find the local maxima of our loaded bedgraph.
void outputLocalMaximalBedGraph(string pathname){

 ifstream inputStream;

    inputStream.open(pathname.c_str());

    if (inputStream.bad()){
        cerr << "WTF??. Error loading in filterBedGraph()";
        abort();
    }
    string lineString;
    stringstream Infostream;
    string start,chr, end;
    float curValue;
    int curPosition, lastPosition;;
    int curEnd, lastEnd;
    float lastValue;
    bool rising=false;

    //Strip and validate header file.
    std::getline(inputStream, lineString);
    Infostream.str(lineString);
    Infostream >>start;

    if (start.find("type=bedGraph")==string::npos)
        {
            cerr << "Failling: No valid Bedgraph header"<<endl;
            abort();
        }
    //Write our new header;
    cout <<"type=bedGraph" <<endl;

    getline(inputStream, lineString);
    Infostream.clear();
    Infostream >>chr;
    Infostream >>lastPosition;
    Infostream >>lastEnd;
    Infostream >>lastValue;

    while(inputStream.eof()!=true){

        getline(inputStream, lineString);
        Infostream.clear();
        Infostream.str(lineString);
        Infostream >> chr;
        Infostream >> curPosition;
        Infostream >> curEnd;
        Infostream >> curValue;
        //Not adjacent, reset
        if (curPosition==(lastPosition+1)){
            //Last position was smaller, we are looking for maxima so continue
            if (lastValue<curValue){
                rising=true;
                lastValue=curValue;
                lastEnd= curEnd;
                lastPosition=curPosition;
            }
            else{
            //Last was a maximal, output
            if (rising){
                    cout << chr << "\t" << lastPosition << "\t" <<  lastEnd<< "\t" << lastValue << endl;
                }
             lastValue=curValue;
             lastEnd= curEnd;
             lastPosition=curPosition;
             rising=false;
            }
        }
        else
        {
            if (rising){
                    cout << chr << "\t" << lastPosition << "\t" <<  lastEnd<< "\t" << lastValue << endl;
                    lastValue=curValue;
                    lastEnd= curEnd;
                    lastPosition=curPosition;
                    rising=false;
                }
           lastPosition= curPosition;
           lastEnd = curEnd;
           lastValue= 0.0f;
        }
    }
}



