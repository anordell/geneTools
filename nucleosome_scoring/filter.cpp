
#include "filter.h"

using namespace std;
void filter(int argc, char* argv[]){
  if(argc <=3)
        {
            cerr<<"Program signature is -b <bedGraph> -s [score threshold=0] -t [top threshold=2] " <<endl;
            return;
        }
        uNucleoBin ourBin;
        float threshold, topthreshold;
        threshold =0;
        topthreshold=2;
        string pathname;
        for (int i = 1; i < argc; i++) {
                if (strcmp(argv[i], "-b")==0) {
                    // We know the next argument *should* be the filename:
                    pathname = argv[i + 1];
                } else if (strcmp(argv[i],"-s")==0) {
                    threshold = atof(argv[i + 1]);
                } else if (strcmp(argv[i],"-t")==0) {
                  topthreshold = atof(argv[i + 1]);
                }
        }

    if ((pathname.size()==0)){
        cerr<<"Program signature is -b <bedGraph> -s [score threshold=0] -t [top threshold=2]" <<endl;
        return;
    }
        filterBedGraph(pathname, threshold,topthreshold);
}


//Load the file and output it filtered. Make generic
//TODO Improve
void filterBedGraph(string pathname, float threshold,float thresholdTop){

    string lineString;
    stringstream Infostream;
    vector<bedScores>::iterator bedit;
    vector<bedScores> bedScoreVec;

    //Strip and validate header file.
    bedScoreVec= loadbedGraph(pathname);

     cout <<"type=bedGraph" <<endl;
    for(bedit=bedScoreVec.begin();bedit!=bedScoreVec.end();bedit++){
        if ((bedit->score>threshold)&&(bedit->score<thresholdTop)){
            cout << bedit->chr<< "\t"<< bedit->position << "\t" << bedit->end << "\t" << bedit->score <<endl;
        }
    }
}
