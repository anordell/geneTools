#include "score.h"


using namespace std;

void score(int argc, char* argv[]){
      if(argc <=3)
        {
            cerr<<"Program signature is -f <filepath> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsComplete=0]";
            return;
        }
        int threshold=0;
        bool isComplete=false;
        string pathname, chromName;// = argv[2];
    for (int i = 1; i < argc; i++) {
                if (strcmp(argv[i], "-f")==0) {
                    // We know the next argument *should* be the filename:
                    pathname = argv[i + 1];
                } else if (strcmp(argv[i],"-c")==0) {
                    chromName = argv[i + 1];
                } else if (strcmp(argv[i],"-n")==0) {
                  threshold = atoi(argv[i + 1]);
                }
                else if(strcmp(argv[i],"-p")==0){
                  int temp = atoi(argv[i+1]);
                    if (temp)
                        isComplete=true;
                }
        }
        //If we did not input file or chrom name
        if ((pathname.size()==0)||(chromName.size()==0)){
            cerr<<"Program signature is -f <filepath> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsPe=0]";
            return;
        }
        uNucleoBin ourBin;
        ourBin.setChromName(chromName);
        cerr << "Trying to Score out file" << endl;
        mainScoring(pathname, ourBin, threshold,isComplete);
        cerr <<"Finished scoring"<< endl;
        ourBin.writebedGraph(cout);
}



void mainScoring(string pathname, uNucleoBin &ourBin, int threshold, bool isComplete){

    uTagsExperiment ourExp;
    uTagsChrom* pChrom;
    ifstream inputStream;

    inputStream.open(pathname.c_str());

    if (inputStream.good()){
        ourExp.loadSamHeader(inputStream);
        }
    else{
        cerr << "WTF?? Loading in mainScoring";
        abort();
    }

    pChrom=ourExp.getpChrom(ourBin.getChromName());


    cerr <<"Treating chromosome of size " << pChrom->getChromSize() << endl;

   // if (isComplete)
   //     ourBin.generateDensityBinFromComplete(&ourExp,pChrom->getChromSize(),threshold);
   // else
    ourBin.generateDensityBin(&ourExp,pChrom->getChromSize(),threshold,isComplete);
    ourBin.generateSMap();
}

