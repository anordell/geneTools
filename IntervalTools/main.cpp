#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include "utility.h"
#include <time.h>
#include "extend.h"

using namespace std;

void printHelp();

int main(int argc, char **argv) {

      string firstArg="";
 if (argc>1)
        firstArg=argv[1];
    if ( firstArg=="-h")
    {
            printHelp();
    }
    else
    if ( firstArg=="extend"){
        extend(argc, argv);
    }
    else
    if (firstArg=="trim"){
       // getDensity(argc, argv);
    }
    else
    if (firstArg=="merge"){
       // parsePhasogram(argc, argv);
    }
     //if (firstArg=="subSetSamFromInterval"){
     //   subSetSamFromInterval(argc, argv);
    }
}

void printHelp(){

  cerr << "Possible argument footprints are:" <<endl;

}
