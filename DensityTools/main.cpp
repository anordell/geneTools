#include "uFormatBase.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include "divideTabFile.h"
#include "getDensity.h"
#include "phasogram.h"
#include "densityFromFile.h"
#include "signalDiff.h"
#include "generateSignal.h"

using namespace std;

void printHelp();

#define SLEEP_LGTH 2  // sleep time in seconds
#define NPOINTS    50 // length of array

int main(int argc, char* argv[])
{

    string firstArg="";
    if (argc>1)
        firstArg=argv[1];

    if ( firstArg=="-h")
    {
        printHelp();
    }
    else if ( firstArg=="divide")
    {
        divideTabFile(argc, argv);
    }
    else if (firstArg=="density")
    {
        getDensity(argc, argv);
    }
    else if (firstArg=="phasogram")
    {
        parsePhasogram(argc, argv);
    }
    else if (firstArg=="densityFromSam")
        densityFromFile(argc, argv);
    else if (firstArg=="signalDiff")
        signalDiff(argc, argv);
    else if (firstArg=="generateSignal")
        generateSignal(argc, argv);
    else
    {
        printHelp();
    }

    return 0;
}

void printHelp()
{

    cerr << "Possible argument footprints are:" <<endl;
    cerr << "divide -f <filepath> -s <Bin Size> -t [default=STRICT; IGNORE/STRICT/EXTEND/ADD]" << endl;
    cerr << "density -f <BinPath> -s <SamPath>" << endl;
    cerr << "phasogram -f <filepath> -s <Graph Size> -p <Pile Size>" << endl;
    cerr << "densityFromSam -s <SamFile> -o [OutputPath]" << endl;
    cerr << "generateSignal " << endl;
}
