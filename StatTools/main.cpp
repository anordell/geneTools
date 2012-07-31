#include "uFormats.h"
#include "uTags.h"
#include "uRegion.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include <time.h>
#include <tclap/CmdLine.h>
#include "peakStats.h"
#include "densityStats.h"
#include "signalMap.h"
#include "diffPlot.h"

void printHelp();
void peakStatistics(int argc, char* argv[]);

using namespace std;
using namespace TCLAP;
int main(int argc, char* argv[])
{

    string firstArg="";
    if (argc>1)
        firstArg=argv[1];

    if ( firstArg=="-h")
    {
        printHelp();
    }
    else if (firstArg=="peakStats")
        peakStatistics(argc, argv);
    else if (firstArg=="densityStats")
        densityStats(argc, argv);
    else if (firstArg=="signalMap")
        signalMap(argc,argv);
    else if (firstArg=="diffPlot")
        diffPlot(argc,argv);
    else
    {
        cerr <<"Invalid arguments. -h for help" <<endl;;
    }
    return 0;

}

void printHelp()
{
    cerr << "Possible argument footprints are: " <<endl;
    cerr << "peakStats " <<endl;
    cerr << "densityStats " <<endl;
}
