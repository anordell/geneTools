#include "NScore.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "uRegion.h"
#include "functions.h"

using namespace std;
struct statsStruct;
void NScore(int argc, char* argv[]){

    uNucleoBin ourBin;


    if(argc <=3)
    {
        cerr<<"Program signature is NScore -f <filepath> -c <chromosomename> -p [IsComplete=0]";
        return;
    }
    int threshold=0;
    bool isComplete=false;
    string pathname, chromName,outputPath;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-f")==0)
        {
            // We know the next argument *should* be the filename:
            pathname = argv[i + 1];
        }
        else if (strcmp(argv[i],"-c")==0)
        {
            chromName = argv[i + 1];
        }
        else if (strcmp(argv[i],"-n")==0)
        {
            threshold = atoi(argv[i + 1]);
        }
         else if (strcmp(argv[i],"-o")==0)
        {
            outputPath = argv[i + 1];
        }
        else if(strcmp(argv[i],"-p")==0)
        {
            int temp = atoi(argv[i+1]);
            if (temp)
                isComplete=true;
        }
    }
    //If we did not input file or chrom name
    if ((pathname.size()==0)||(chromName.size()==0))
    {
        cerr<<"Program signature is Nscore -f <filepath> -c <chromosomename>  -p [IsPe=0]";
        return;
    }

    //vector<uRegion> ourRegions;
    //Load our Sam File
    {
        uTagsExperiment ourExp;
        ifstream inputStream;
        inputStream.open(pathname.c_str());
        //Load our Sam file
        if (inputStream.good())
            ourExp.loadSamHeader(inputStream);
        else
        {
            cerr << "Error loading in NScore";
            abort();
        }
        //Loaded, generate our Density map per bp
        ourBin= loadDensityFromTags(ourExp, threshold,chromName, isComplete);
        cerr << "Finished density bins, Starting map Regions" <<endl;
        //Make our first map based on Nature score
        ourBin.generateSMap();
        statsStruct ourStats=ourBin.getSD();

        /**< Secondary options */
     //   auto ourRegions = makeRegionsfromBin(ourBin,ourExp);
        //Get SD and Mean
     //   statsStruct ourSD=returnSdfromRegions(ourRegions);
        /**< Erase */

        cerr <<"Measuring Sum and Mean" << endl;
        cerr <<"SD is" <<ourStats.sd<< endl;
        cerr <<"Mean is" <<ourStats.mean<< endl;

       // utility::pause_input();


        if (outputPath.size()!=0)
        {
            ofstream outputOS(outputPath);
            ourBin.writebedGraph(outputOS);
        }


        ourBin.GaussianNormSMap(ourStats.sd, ourStats.mean);

        ourBin.writebedGraph(cout);





    }


}
