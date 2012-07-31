#include "uFormats.h"
#include "uTags.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include "utility.h"
#include <time.h>

using namespace std;

int main(int argc, char **argv) {



    string pathname;


    pathname = argv[1];
    uTagsExperiment ourExp;
    uTags tempTag;
    ifstream inputStream;

    inputStream.open(pathname.c_str());
    if (inputStream.good())
        ourExp.loadFromSam(inputStream);
    else{
        cerr << "WTF?? Loading";
        abort();
    }

    ofstream plusStream("fTags.txt");
    ofstream minusStream("rTags.txt");

     ourExp.applyOnSites([&](uTags & Elem)
    {
    if (Elem.getChr()=="chr21"){
       if ( Elem.getStrand() =='+'){
            plusStream<< Elem.getStart()<<endl;
        }
        else
        {
            minusStream<<Elem.getEnd()<<endl;
        }
        }
    }

    );




}
