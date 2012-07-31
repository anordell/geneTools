
#include "PEComplete.h"


using namespace std;
void completePeSam(string pathname, ostream& output){

    uTagsExperiment ourExp;
    uTags tempTag;
    ifstream inputStream;

    inputStream.open(pathname.c_str());

    if (inputStream.good())
        //ourExp.loadSamHeader(inputStream);
        ourExp.loadFromSam(inputStream);
    else{
        cerr << "WTF?? Loading";
        abort();
    }

    ourExp.writeCompletedPESamToOutput(output);

}


