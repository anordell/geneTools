
#include "uTags.h"
#include "SamTrim.h"
using namespace std;
void SamTrim(int trim, string pathname, ostream& output){

    SamTrim(trim, trim,pathname, output);
}



void SamTrim(int left, int right, string pathname, ostream& output){

 uTagsExperiment ourExp;
    uTags tempTag;
    ifstream inputStream;

    inputStream.open(pathname.c_str());

    if (inputStream.good())
        ourExp.loadFromSam(inputStream);
        //ourExp.loadFromSam(inputStream);
    else{
        cerr << "WTF?? Loading";
        abort();
    }

    ourExp.writeTrimmedSamToOutput(output, left, right);


}
