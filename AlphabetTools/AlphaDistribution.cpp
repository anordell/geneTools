#include <tclap/CmdLine.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include <time.h>
#include <random>
#include <string.h>
#include "utility.h"
#include "uAlphabet.h"
#include "AlphaDistribution.h"



using namespace std;
const char HEADERSTART='>';
void AlphaDistribution(int argc, char* argv[])
{

    try
    {
        TCLAP::CmdLine cmd("Analyse of an fa file, spitting out alphabet distribution per sequence and globaly", ' ', "0.1");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : AlphaDistribution", false," ","String");
        TCLAP::ValueArg<std::string> fileArg("f","file","fa file containing our data ",true," ","string", cmd);
        TCLAP::MultiArg<char> IDArg("I","Ignore","Letters to ignore for local distributions",false,"char",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Path to output file",false,"","string",cmd);


        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */
        // auto filepath = fileArg.getValue();
        auto IgnoreCharVector = IDArg.getValue();



        ifstream inputStream;
        utility::loadStream(fileArg.getValue(), inputStream);

        auto alphaData = loadData(inputStream);
        processMap(alphaData,IgnoreCharVector);

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (string a)
    {
        cerr <<"Some kind of error string " << a <<endl;

    }
    catch (...)
    {
        cerr <<"Some kind of error" <<endl;

    }


}
/**< Add checks to this */
AlphaMap loadData(ifstream & inputStream)

{
    cerr << "Loading" <<endl;
    AlphaMap returnMap;

    string templine;
    string curHeader;
    /**< Parse our file, header then data */
    while(!std::getline(inputStream,templine).eof())
    {
        if (templine.at(0)==HEADERSTART)
        {
            cerr << "Writing header as" << templine <<endl;
            curHeader=templine;
        }
        else
for(char curLetter:templine)
            {
                returnMap[curHeader][curLetter]++;
            }
    }
    return returnMap;
}

void processMap(const AlphaMap & alphaData, vector<char> ignoreVecChar)
{
    stringstream outputString;
    float globalCount=0, globalCountNoIgnore=0;
    map<char,int> globalMap;
    /**< For every sequence */

    cerr << "Writing alphaData of size "<< alphaData.size() <<endl;

    for (auto seqMap: alphaData)
    {
        int localCount=0, localCountNoIgnore=0;
        /**< Every char of that sequence */
        for(auto charPair:seqMap.second )
        {
            /**< Add to global count */
            globalCount+=charPair.second;
            if ( find(ignoreVecChar.begin(),ignoreVecChar.end(),charPair.first )==ignoreVecChar.end())
            {
               // cerr << "found"<<endl;
                globalCountNoIgnore+=charPair.second;
            }

            localCount+=charPair.second;
            if ( find(ignoreVecChar.begin(),ignoreVecChar.end(),charPair.first )==ignoreVecChar.end())
                localCountNoIgnore+=charPair.second;

            globalMap[charPair.first]+=charPair.second;

        }
        /**< Output local data for that sequence */
        outputString << seqMap.first <<endl;
        outputString<< localCountNoIgnore<< " bp / " <<localCount << endl;
        outputString<< "Percent covered by all  = " << (float)localCountNoIgnore/(float)localCount << endl;
        for(auto charPair:seqMap.second )
        {
            outputString<<charPair.first <<endl;
            outputString<<charPair.second<< " bp / " <<localCountNoIgnore << endl;
            outputString<<(float)charPair.second/(float)localCountNoIgnore << endl;
        }
    }
    ofstream outputStream("test.txt");
    outputStream << outputString.str();
}
