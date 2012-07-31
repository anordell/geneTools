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
#include "utility.h"
#include <time.h>
#include <sstream>
using namespace std;

typedef std::tuple <std::string, int, int> chr_tuple;
void generateRandomExpFromBed(int argc, char* argv[]);
void generateRandomExp(int argc, char* argv[]);
vector<chr_tuple> readChrSizes(ifstream& input);
std::tuple <map<string, int>, float, float> readBed(ifstream& input, const vector<chr_tuple> ourVec );
void writeHelp();
int main(int argc, char* argv[])
{

    string chr;
    bool fromBed=false, help=false;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-chrFile")==0)
        {
            // We know the next argument *should* be the filename:
            fromBed=true;
        }
        else if (strcmp(argv[i],"-bedFile")==0)
        {
            fromBed=true;
        }
        else if (strcmp(argv[i],"-help")==0)
        {
            help=true;
        }
    }

    if ((help)||(argc==1))
    {
        writeHelp();
        return 0;
    }
    if (fromBed)
        generateRandomExpFromBed(argc, argv);
    else
        generateRandomExp(argc, argv);

}

void writeHelp()
{

    cerr<<"Program signature is -chrsize <chromsize> -chr <chromosomename> -count <samplesize> -regsize <size of random regions> opt: -exclusionFile [exclusion file] -seed [seednumber] -id [id name for batch] -sd [standard deviation of size]" << endl;
    cerr << " or " << endl;
    cerr<< "Program signature is -bedFile <bed File Path> -chrFile <Chromosome Info Path>  opt: -seed [seed] -id [Id for batch]"<<endl;

    return;
}


void generateRandomExp(int argc, char* argv[])
{

    string chr,ID, exclPath;
    int seed=0,count=0, regionsize=0, sd=0;
    int chrsize2=0;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-chrsize")==0)
        {
            chrsize2 = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"-chr")==0)
        {
            chr = argv[i + 1];
        }
        else if (strcmp(argv[i],"-seed")==0)
        {
            seed = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"-regsize")==0)
        {
            regionsize = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"-count")==0)
        {
            count = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"-sd")==0)
        {
            sd = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"-id")==0)
        {
            ID = argv[i+1];
        }
        else if (strcmp(argv[i],"-exclusionFile")==0)
        {
            exclPath= argv[i+1];
        }
    }
    if ((chrsize2==0)||(chr.size()==0)||(count==0)||(regionsize==0))
    {
        cerr<<"Program signature is -chrsize <chromsize> -chr <chromosomename> -count <samplesize> -regsize <size of random regions> opt: -exclusionFile [exclusion file] -seed [seednumber] -id [id name for batch] -sd [standard deviation of size]";
        return;
    }

    if (chrsize2<=regionsize)
    {
        cerr<<"Input -chrsize larger then -regsize";
        return;
    }

    if (ID.size()==0)
    {
        cerr<<"No id specified, defaulting" << endl;
        ID="rnd";
    }

    string inputname, outputname;
    // int start, stop;
    ofstream outputfile;
    ifstream inputfile;

    uTagsChrom newTags;

    if (seed==0)
    {
        // seed=srand(time(NULL));
        //Not a good random, but hey
        seed = (time(NULL)%10000);
        cerr <<"No seed or seed=0 specified, using " << seed << endl;
    }

    uTagsChrom exclExp;
    uRegionChrom exclReg;
    std::stringstream sstm;
    if (exclPath.size())
    {
        ifstream exInput(exclPath.c_str());
        if (exInput.fail())
        {
            cerr<<"Error opening file "<< exclPath.c_str()<<endl;
            return;
        }

        std::string exChr, buffer;
        int exStart,exEnd;
        while(!exInput.eof())
        {
            std::getline(exInput, buffer);
            sstm << buffer;
            sstm >> exChr;
            sstm >> exStart;
            sstm >> exEnd;
            sstm.clear();
           // cerr << "Zone" <<"\t"<< exChr << "\t" << exStart << "\t" << exEnd << endl;
            // std::getline(exInput, exChr, '\t').eof();
            uTags excReg(exChr,exStart,exEnd);

            if (exChr==chr)
                exclExp.addSite(excReg);
        }


    }
    newTags.setChr(chr);
    newTags.setChromSize(chrsize2);

    ofstream haha;
//    haha.open("ourExclusionlist.txt");
//    exclExp.outputBedFormat(haha);
    exclExp.sortSites();

  //  exclExp.outputBedFormat(cout);
    while(newTags.count()!=count)
    {
        std::mt19937 engine;
        //seed+=5;
        engine.seed(seed);
        uTagsChrom temptags;

        newTags.addNRandomSite(regionsize, count,engine,exclExp, sd, ID);
    }

    newTags.outputBedFormat(cout);

    cerr << "Seed used was : "<< seed<< " . Please note this if you wish to reproduced results";

    return;
}


void generateRandomExpFromBed(int argc, char* argv[])
{

    vector<chr_tuple> chrSizes;
    string chrFilePath, bedFilePath, ID;
    int seed=0;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i],"-chrFile")==0)
        {
            chrFilePath = argv[i + 1];
        }
        else if (strcmp(argv[i],"-seed")==0)
        {
            seed = atoi(argv[i+1]);
        }
        else if (strcmp(argv[i],"-bedFile")==0)
        {
            bedFilePath = argv[i+1];
        }
        else if (strcmp(argv[i],"-id")==0)
        {
            ID = argv[i+1];
        }
    }
    if ((bedFilePath.size()==0)||(chrFilePath.size()==0))
    {
        cerr<< "Program signature is -bedFile <bed File Path> -chrFile <Chromosome Info Path>  opt: -seed [seed] -id [Id for batch]"<<endl;
        // cerr<<"If -bedFile or -chrFile is used, both must be used.";
        return;
    }

    cerr<<"-bedFile & -chrFile specified, ignoring other parameters";
    if (ID.size()==0)
    {
        cerr<<"No id specified, defaulting" << endl;
        ID="rnd";
    }
    if (seed==0)
    {
        //   seed=srand(time(NULL));
        //Not a good random, but hey
        seed = (time(NULL)%100000);
        cerr <<"No seed or seed=0 specified, using " << seed << endl;
    }

    ifstream bedFile, chrFile;
    bedFile.open(bedFilePath.c_str());
    chrFile.open(chrFilePath.c_str());
    if ((bedFile.bad())||(chrFile.bad()))
    {
        cerr<<"Error opening file from either "<< bedFilePath << " or " << chrFilePath<<endl;
        return;
    }

    //  cerr << bedFilePath ;
    std::mt19937 engine;
    engine.seed(seed);
    chrSizes=readChrSizes(chrFile);

    std::tuple <map<std::string, int>, float, float> mapTuple = readBed(bedFile,chrSizes);
    float mean    = get<1>(mapTuple);
    float sd      = get<2>(mapTuple);

    std::map<std::string, int> chrMap = get<0>(mapTuple);
    string curchr;
    chr_tuple tempChr;
    int chrsize, entrycount=0;
    //for(auto it=chrSizes.begin(); it != chrSizes.end(); it++ ){
for (auto& it : chrSizes)
    {
        tempChr= it;
        curchr = get<0>(tempChr);
        chrsize = get<1>(tempChr);
        {
            uTagsChrom newTags;
            //Pour chaque chromosome, vérifier si il faut en créer
            if(chrMap[curchr]>0)
            {
                newTags.setChr(curchr);
                newTags.setChromSize(chrsize);
                cerr << mean << " " << sd << endl;

                while(newTags.count()!=chrMap[curchr])
                {
                    newTags.addNRandomSite(mean, chrMap[curchr],engine, sd,ID);

                }
                newTags.outputBedFormat(cout);

                entrycount+=mean;
            }
        }
    }
    cerr << "Generated " <<entrycount << " lines with a mean of " << mean << " and sd of " <<sd << endl;
    cerr << "Seed used was : "<< seed<< " . Please note this if you wish to reproduced results";

    return;
}

// Get our chromosome sizes from our bed file
// Return them in our brand new tupples.
vector<chr_tuple> readChrSizes(ifstream& input)
{
    string tempChr;
    int tempsize;
    chr_tuple tempTuple;
    vector<chr_tuple> returnValues;

    while (input.eof()==false)
    {

        input>>tempChr;
        input>>tempsize;
        get<0>(tempTuple)=tempChr;
        get<1>(tempTuple)=tempsize;
        returnValues.push_back(tempTuple);

        input.peek();
        if (input.eof())
            break;
    }

    return returnValues;
}

// Get our count per chromosome
// Return them as a map
//and return SD
std::tuple <map<std::string, int>, float, float> readBed(ifstream& input, const vector<chr_tuple> ourVec )
{

    string tempChr,string, chr,temp;
    int start, end, count=0, totsize=0;
    chr_tuple tempTuple;
    vector <float> sizes;
    std::map<std::string, int> ourMapCount;
//vector<chr_tuple> returnValues;
    std::string curLine;

    while (!std::getline(input, curLine).eof())
    {
        {
            utility::Tokenizer data(curLine);
            //   cerr <<curLine<< endl;
            data.NextToken();
            tempChr = data.GetToken();
            // cerr << tempChr << endl;
            data.NextToken();
            start   = utility::stringToInt(data.GetToken());
            data.NextToken();
            end     = utility::stringToInt(data.GetToken());;

            count++;
            ourMapCount[tempChr]=(ourMapCount[tempChr]+1);
            totsize += (end-start);
            sizes.push_back(start-end);
        }
        //input.peek();
        if (input.eof())
            break;
    }
    float mean = totsize/count;
    float sd = utility::getSd(sizes, mean);


    return tie(ourMapCount,mean,sd);
//return std::tuple <map<std::string, int>, float, float> (ourMapCount, mean, sd);
}
