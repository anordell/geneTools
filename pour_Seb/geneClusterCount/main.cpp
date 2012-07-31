#include <iostream>
#include "gnuplot_i.hpp"
#include "include/comparisonGroup.h"
#include "utility.h"
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>

using namespace std;


vector<string> parseHeader(std::string pHeader);
void exploreTree(comparisonGroup& currentGroup, std:: vector<comparisonGroup>::iterator itrCurBasic, std:: vector<comparisonGroup>::iterator itrCurEnd, vector<comparisonGroup> & resultGroup);
int main(int argc, char* argv[])
{

    //Load our initial Data
    string filename;
    string outputPrefix;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-f")==0)
        {
            filename = argv[i + 1];
        }
        else if (strcmp(argv[i],"-s")==0)
        {
            outputPrefix = argv[i + 1];
        }
    }

    if ((filename.size()==0)||(outputPrefix.size()==0))
    {
        cerr<<"Program signature is -f <filepath> -s <Output Prefix>";
        return 0;
    }

    ifstream inputStream(filename);

    if (!(inputStream.good()))
    {
        cerr << "Unable to load. Vincent is to blame!";
        abort();
    }

    string header,line;
    getline (inputStream,header);

    vector<string> geneList = parseHeader(header);
    vector<comparisonGroup> initialExperiment;
    initialExperiment.resize(geneList.size());
    for (unsigned int i=0; i< geneList.size(); i++)
    {
        initialExperiment.at(i).addGeneID(geneList.at(i));
    }

    //Create one Cluster per gene and load corresponding data.
    while ( !std::getline(inputStream, line).eof() )
    {

        {
            utility::Tokenizer lineParser(line);
            lineParser.NextToken();
            string ClusterID=lineParser.GetToken();
            int column=0;
            while(lineParser.NextToken())
            {
                if (lineParser.GetToken()=="1")
                    initialExperiment.at(column).addClusterData(ClusterID,PRES);
                else
                    initialExperiment.at(column).addClusterData(ClusterID,ABS);
                column++;
            }
        }
    }

for ( comparisonGroup curComparison : initialExperiment )
    {
        cout << curComparison.getConcatGeneID()<< " " << curComparison.getPresentCount() << " " << curComparison.getAbsentCount() << " "  <<curComparison.getMismatchCount() << endl;
    }



    vector<comparisonGroup> resultGroup;

    std:: vector<comparisonGroup>::iterator itrCur;



    auto itrEnd = initialExperiment.end();
    for ( auto itr = (initialExperiment.begin()); itr!=initialExperiment.end(); itr++)
    {

        exploreTree(*itr, itr,itrEnd, resultGroup);
    }


for ( comparisonGroup curComparison : resultGroup )
    {
        cout << curComparison.getGeneCount()<< " " << curComparison.getPresentCount() << " " << curComparison.getAbsentCount() << " "  <<curComparison.getMismatchCount() << endl;
    }


    return 0;
};

void exploreTree(comparisonGroup& currentGroup, std:: vector<comparisonGroup>::iterator itrCurBasic, std:: vector<comparisonGroup>::iterator itrCurEnd, vector<comparisonGroup> & resultGroup)
{

    if (itrCurBasic!=itrCurEnd)
        itrCurBasic++;
    while(itrCurBasic!=itrCurEnd)
    {
        //new element
        comparisonGroup grtoAdd;
        grtoAdd=currentGroup.returnCompResult(*itrCurBasic);
        resultGroup.push_back(grtoAdd);
        exploreTree(grtoAdd,itrCurBasic,itrCurEnd,resultGroup);
        itrCurBasic++;
    }
}


vector<string> parseHeader(std::string pHeader)
{
    using namespace utility;

    utility::Tokenizer tokenHeader(pHeader);
    vector<string> returnVector;
    //Discard header
    tokenHeader.NextToken();

    //Every gene name into the return vector
    while(tokenHeader.NextToken())
        returnVector.push_back(tokenHeader.GetToken());


    return returnVector;
}

