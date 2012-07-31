
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <random>
#include <string.h>
#include "utility.h"
#include "generateEpiAlpha.h"
#include "AlphaDistribution.h"

using namespace std;

void printHelp();

int main(int argc, char* argv[])
{

    string firstArg="";
    if (argc>1)
        firstArg=argv[1];

    if ( firstArg=="-h")
    {
        printHelp();
    }
    else if (firstArg=="generateEpiAlpha")
        generateEpiAlpha(argc, argv);
    else if (firstArg=="AlphaDistribution")
        AlphaDistribution(argc, argv);
    else
    {
        cerr <<"Invalid arguments. -h for help" <<endl;;
    }

    return 0;

}


void printHelp()
{

    cerr << "Possible tools calls are: " <<endl;
    cerr << "generateEpiAlpha  " <<endl;
    cerr << "AlphaDistribution" <<endl;
    cerr << "<tool> -h for help  " <<endl;
}

