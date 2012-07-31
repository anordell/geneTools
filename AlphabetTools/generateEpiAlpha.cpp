#include <tclap/CmdLine.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include <random>
#include <string.h>
#include "utility.h"
#include "uAlphabet.h"
#include "generateEpiAlpha.h"


using namespace std;

const string tab="\t";


bool allChrom(uAlphabetChrom value)
{
    if (value.count()>0)
        return true;
    else
        return false;
}


std::function<bool(uAlphabetChrom)> chromPred = &allChrom;


void generateEpiAlpha(int argc, char* argv[])
{

    try
    {
        TCLAP::CmdLine cmd("Create an epigenetic alpabet based on N marcs", ' ', "0.3");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : generateEpiAlpha", false," ","String");
        TCLAP::MultiArg<std::string> fileArg("f","file","Tab file containing Epigenetic region",true,"string",cmd);
        TCLAP::MultiArg<std::string> IDArg("I","ID","Identifier for each Epigenetic region",true,"string",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Path to output file",false,"","string",cmd);
        //  TCLAP::SwitchArg isBedSwitch("b","bed","If a bed file", cmd, false);
        //     TCLAP::SwitchArg isSamSwitch("s","sam","If a Sam file", cmd, false);
        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */
        auto fileVector = fileArg.getValue();
        auto IDVector = IDArg.getValue();

        /**< Need same number of ID and files */
        if (fileVector.size()!=IDVector.size())
            throw 10;

        /**< Const genome size for Hg18 */
        map<string, int> sizeMap;
        sizeMap["chr1"]=247249719;
        sizeMap["chr2"]=242951149;
        sizeMap["chr3"]=199501827;
        sizeMap["chr4"]=191273063;
        sizeMap["chr5"]=180857866;
        sizeMap["chr6"]=170899992;
        sizeMap["chr7"]=158821424;
        sizeMap["chr8"]=146274826;
        sizeMap["chr9"]=140273252;
        sizeMap["chr10"]=135374737;
        sizeMap["chr11"]=134452384;
        sizeMap["chr12"]=132349534;
        sizeMap["chr13"]=114142980;
        sizeMap["chr14"]=106368585;
        sizeMap["chr15"]=100338915;
        sizeMap["chr16"]=88827254;
        sizeMap["chr17"]=78774742;
        sizeMap["chr18"]=76117153;
        sizeMap["chr19"]=63811651;
        sizeMap["chr20"]=62435964;
        sizeMap["chr21"]=46944323;
        sizeMap["chr22"]=49691432;
        sizeMap["chrX"]=154913754;
        sizeMap["chrY"]=57772954;



        std::string outputPath = outputArg.getValue();
        /**< Name treatment */
        vector<uAlphabetExperiment> alphabetData;

        /**< Load our genome data */
        alphabetData=loadData(fileVector, IDVector, sizeMap);

        vector<int> IDVec;
        /**< Set and generate ID list */
        IDVec = setID(alphabetData);

        /**< Debug */
        cerr << "ID List is : ";
        for(int ID: IDVec)
            cerr <<ID <<tab;
        cerr <<endl;


        auto IDMap = createAlphabet(IDVec);
        IDMap[0]='a';


                ofstream index("Index.txt");
                /**< Debug */
                index << "Alphabet is : ";
                for(auto ID: IDMap){

                    index <<ID.first <<tab<<ID.second <<endl;

                }


                index << "Basic Index is : "<<endl;
                for(uAlphabetExperiment curAlph:alphabetData)
                {
                   index << curAlph.getID() << tab << curAlph.getExplain() <<endl;
                }



        auto genomeMap =overlapGenomeWithMarks(alphabetData,IDMap,outputPath);





    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (int e)
    {

        if (e==10)
            std::cerr << "Please input an equal number of -f and -I parameters"<< std::endl;

        if (e==20)
            std::cerr << "Failed loading data"<< std::endl;
    }



}

/** \brief Write our calculated Genomic alphabet
 *
 * \param map<string, vector<int>> genomeMap Map with our genome Data
 * \param out ofstream& Where to write it.
 * \return void
 *
 */
void writeMap(map<string,vector<int>> genomeMap, std::ostream& out)
{

for(auto & mapItem:genomeMap )
    {
        out<<    mapItem.first;
for (auto vecItem : mapItem.second)
            out<< tab<<   vecItem ;
        out << endl;
    }


}

/** \brief Load our data
 *
 * \param pathVector vector<string> Vector containing path to every genomic mark
 * \param IDVector vector<string> Vector containing the ID/Explanation of every mark
 * \return vector<uAlphabetExperiment> Our Loaded data
 *
 */
vector<uAlphabetExperiment> loadData(vector<string> pathVector, vector<string> IDVector,map<string,int> sizeMap)
{
    vector<uAlphabetExperiment> resultVector;

    /**< Load data and set size */
    try
    {
        int i=0;
        for(string path:pathVector)
        {
            {
                ifstream stream;
                uAlphabetExperiment currentExp;
                utility::loadStream(path,stream);
                currentExp.loadFromTabFile(stream);
                currentExp.setExplain(IDVector.at(i));
                for(auto chromSize:sizeMap)
                    currentExp.setChrSize(chromSize.first, chromSize.second);

                resultVector.push_back(currentExp);
            }
            i++;
        }
    }
    catch(string a)
    {
        cerr <<"Failed loading file in loadData " <<  a<<endl;
        throw 20;
    }



    return resultVector;

}

/** \brief Add each mark to to genome. Sum the value of every mark, 0 being no epigenetic mar.
 *
 * \param p_markData vector<uAlphabetExperiment>
 * \return map<string,vector<int>>
 *
 */
map<string,vector<int>> overlapGenomeWithMarks(vector<uAlphabetExperiment> p_markData, std::map<int,char> AlphabetMap, string outputPath)
{


    uAlphabetExperiment* markdata=&(p_markData.at(0));
    map<string,vector<int>> returnMap;

     /**< Erase output file if any */
    {
     if (outputPath.size()!=0)
         ofstream outputOS(outputPath);
    }


    for (auto pchrom = markdata->first(); pchrom!=markdata->last(); pchrom++)
    {
        {
            vector<int> curVector;
            curVector.resize(pchrom->second.getChromSize());
            cerr <<pchrom->first << endl;
            cerr <<curVector.size() << endl;

            for(auto& Experiment: p_markData )
            {
                Experiment.modifyAlphabet(pchrom->first,curVector);
            }

            //curVector
            if (outputPath.size()!=0)
            {

                ofstream outputOS(outputPath, ios::app);
                 outputOS <<">" <<pchrom->first<<endl;
                    for (auto value: curVector)
                        outputOS<<AlphabetMap[value];
                outputOS<<endl;

                // writeMap(genomeMap, outputOS);
            }
            else
            {
                cout <<">" << pchrom->first<<endl;
                    for (auto value: curVector)
                    cout<<AlphabetMap[value];
                cout<<endl;
            }

            //returnMap.insert(pair<string,vector<int>>(pchrom->second.getChr(),curVector));
        }
    }

    /* for(auto& Experiment: p_markData ){
         Experiment.modifyAlphabet(returnMap);
     }*/

    return returnMap;
}


/** \brief Assign a unique base 2 ID to each mark
 *
 * \param histList vector<uAlphabetExperiment>&  List of loaded marks
 * \return vector<int> List of assigned ID
 *
 */
vector<int> setID(vector<uAlphabetExperiment> & histList)
{
    vector<int> result;
    int cur=1;
for(uAlphabetExperiment& histExp : histList )
    {
        histExp.setID(cur);
        result.push_back(cur);
        cur=cur*2;
    }
    return result;
}


/** \brief From a list of ID, assign a list of characters to each combination
 *
 * \param idList const vector<int>& List of numerical ID
 * \return map<int,char> Map linking each value to a char.
 *
 */
map<int,char> createAlphabet(const vector<int> & idList)
{

    map<int,char> alphabet;
    char currentChar='b';
    auto itrBegin = idList.begin();
    auto itrEnd = idList.end();
    exploreAlphabet(alphabet,itrBegin,itrEnd, currentChar, 0);
    return alphabet;
}

/** \brief Recursive call for createAlphabet. Generate every permuation
 *
 * \param map<int
 * \param alphabet char>&
 * \param itrCurBasic std::vector<int>::const_iterator
 * \param itrCurEnd std::vector<int>::const_iterator
 * \param curChar char&
 * \param curValue int
 * \return void
 *
 */
void exploreAlphabet(map<int,char>& alphabet,  std::vector<int>::const_iterator itrCurBasic,  std::vector<int>::const_iterator itrCurEnd, char & curChar, int curValue)
{
    while(itrCurBasic!=itrCurEnd)
    {
        //new element
        int newcurValue=curValue+(*itrCurBasic);
        alphabet.insert(pair<int,char>(newcurValue,curChar));
        curChar++;
        itrCurBasic++;
        exploreAlphabet(alphabet,itrCurBasic,itrCurEnd,curChar,newcurValue);
    }
}
