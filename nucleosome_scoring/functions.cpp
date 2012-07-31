/*
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "uTags.h"
#include "uRegion.h"
#include "functions.h"

using namespace std;


#define REGSIZE 30

std::vector<bedScores> loadbedGraph(string bedpath)
{
    std::ifstream in(bedpath);

    vector<bedScores> returnVec;
    returnVec= move(loadbedGraph(in));

    return returnVec;
}



std::vector<bedScores> loadbedGraph(std::ifstream& inputStream)
{
   // ifstream inputStream;

   // inputStream.open(bedpath.c_str());
    if (inputStream.bad())
    {
        cerr << "WTF??. Error loading in loadbedGraph()";
        abort();
    }

    string lineString,header;
    stringstream Infostream;
    std::getline(inputStream, lineString);
    Infostream.str(lineString);
    Infostream >>header;

    if (header.find("type=bedGraph")==string::npos)
    {
        cerr << "Error, this is not a valid bedgraph file"<<endl;
        abort();
    }
    bedScores tempscore;
    vector<bedScores> returnVector;
    while(inputStream.eof()!=true)
    {

        getline(inputStream, lineString);
        //Reset error flags our nothing work.
        Infostream.clear();
        Infostream.str(lineString);
        Infostream >>tempscore.chr;
        Infostream >>tempscore.position;
        Infostream >>tempscore.end;
        Infostream >>tempscore.score;
        returnVector.push_back(tempscore);
    }
    return returnVector;
}

//Pass a list of maxima to this function and it will select specific peaks within competing regions
//As we only want to remove one potentional nucleosome per position.
vector<uRegion> filterNucleoMaxima(uTagsExperiment* ourTagExp, string chr, vector<bedScores> &vecScores)
{

    uRegion tempRegion;
    tempRegion.setChr(chr);
    vector<bedScores>::iterator it;
    vector<int>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    vector<uRegion> regionsKept, allregions;
    // load our maxima
    cerr << " We have " << vecScores.size() << endl;
    cerr << "Loading our maxima" << endl;
    for ( it=vecScores.begin(); it < vecScores.end(); it++ )
    {

        tempRegion.setEnd(it->position+30);
        tempRegion.setStart(it->position-30);
        tempRegion.setScore(it->score);
        //TODO fix this
        tempRegion.setCount( ourTagExp->getSubsetCount(tempRegion,OverlapType::OVERLAP_CENTRE));
        allregions.push_back(tempRegion);
    }

    cerr << "Checking our regions" << endl;
    for ( regit=allregions.begin(); regit < allregions.end(); regit++ )
    {
        bool ok=true;
        ic= regit;
        //Go backwards until before the region
        while ((ic->getEnd()>=regit->getStart())&&(ic!=allregions.begin()))
        {
            ic--;
        }
        //Go forward until past
        while ((ic->getStart()<=regit->getEnd())&&(ic!=allregions.end()))
        {
            if (ic->getScore()>regit->getScore())
            {
                ok =false;
                break;

            }
            ic++;
        }

        if (ok)
            regionsKept.push_back(*regit);
    }

    return regionsKept;
}


uTagsExperiment decompPass(uTagsChrom* ourTagChrom, vector<uRegion> vecRegions)
{

    uTagsExperiment returnTagsExp;
    uTagsChrom tempChrom;
    uRegion tempRegion;
    vector<uRegion>::iterator regit;
    vector<uRegion> regionsRemoved;
    //For each maxima, create a region, subset the tags, etc.
    for ( regit=vecRegions.begin(); regit < vecRegions.end(); regit++ )
    {
        tempChrom=ourTagChrom->getSubset(regit->getStart(), regit->getEnd(), OverlapType::OVERLAP_CENTRE);
        returnTagsExp.combine(tempChrom);
    }
    //Remove all our tags from the original dataSet
    int count=0;
    for ( regit=vecRegions.begin(); regit < vecRegions.end(); regit++ )
    {
        count++;
        ourTagChrom->removeSubset(regit->getStart(), regit->getEnd(),OverlapType::OVERLAP_CENTRE);
        //   cerr << "Completed region " << count << endl;
    }
    return returnTagsExp;
}


void decomp(int argc, char* argv[])
{

    if(argc <=5)
    {
        cerr<<"Program signature is -s <samFile> -c <chromosone> -b <BedGraph Maxima> -o <SamLeft> <SamDecomp>" <<endl;
        return;
    }
    string sampath;
    string chromName;
    string bedpath;//; = argv[4];
    string samLeft;//; = argv[5];
    string samDec;// = argv[6];
    uTagsExperiment ourExp;
    uTagsChrom* pChrom;
    ifstream inputStream;
    ifstream bedStream;
    ofstream ofLeft, ofDec;
    stringstream Infostream;
    vector<bedScores> scoreVec;
    bedScores tempScore;

//Parse parameters
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-s")==0)
        {
            sampath = argv[i + 1];
        }
        else if (strcmp(argv[i],"-c")==0)
        {
            chromName = argv[i + 1];
        }
        else if (strcmp(argv[i],"-b")==0)
        {
            bedpath= argv[i + 1];
        }
        else if(strcmp(argv[i],"-o")==0)
        {
            samLeft = argv[i + 1];
            samDec = argv[i + 2];
        }
    }
//Validate we loaded them all
    if ((sampath.size()==0)||(chromName.size()==0)||(bedpath.size()==0)||(samLeft.size()==0)||(samDec.size()==0))
    {
        cerr<<"Program signature is -s <samFile> -c <chromosone> -b <BedGraph Maxima> -o <SamLeft> <SamDecomp>" <<endl;
        return;
    }

    inputStream.open(sampath.c_str());
    if (inputStream.good())
        ourExp.loadFromSam(inputStream);
    //ourExp.loadFromSam(inputStream);
    else
    {
        cerr << "WTF?? Loading";
        abort();
    }
    cerr << "Loading Maxima file" << endl;
    bedStream.open(bedpath.c_str());
    if (bedStream.good())
    {
        string chr;
        string lineString;
        int start, end;
        float score;
        //Skip header
        bedStream>>chr;
        while(bedStream.eof()!=true)
        {
            getline(bedStream, lineString);
            Infostream.clear();
            Infostream.str(lineString);
            Infostream >> chr;
            Infostream >> start;
            Infostream >> end;
            Infostream >> score;

            tempScore.position=start;
            tempScore.score=score;
            scoreVec.push_back(tempScore);
        }
    }
    else
    {
        cerr << "WTF?? Loading on Bed";
        abort();
    }
    //  uTagsChrom* pChrom;
    // vector<int> densityMap;
    pChrom=ourExp.getpChrom(chromName);
    vector<uRegion> regionVector;
    //createDensityMap(pChrom,densityMap);
    cerr << "Filtering our maxima" << endl;
    regionVector=filterNucleoMaxima(&ourExp, chromName,scoreVec);

    ofLeft.open(samLeft.c_str());
    ofDec.open(samDec.c_str());

    uTagsExperiment expDec;
    cerr << "Decomp our data" << endl;
    expDec=decompPass(pChrom,regionVector);
    expDec.sortData();
    ourExp.writeSamToOutput(ofLeft);
    expDec.writeSamToOutput(ofDec);
}

statsStruct returnSdfromRegions(vector<uRegion> regionVec)
{

    vector<bedScores>::iterator it;
    vector<float>::iterator floatit;
    vector<int>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    int sum, mean, sd,sumsq;
    sum=0;
    sumsq=0;
    //Establish our score passed on density average and sigma
    vector<int> deviations;
    for ( regit=regionVec.begin(); regit < regionVec.end(); regit++ )
    {
        sum +=regit->getCount();
    }

    mean = sum/regionVec.size();
    for ( regit=regionVec.begin(); regit < regionVec.end(); regit++ )
    {
        deviations.push_back(regit->getCount()-mean);
    }
    for (intit=deviations.begin(); intit < deviations.end(); intit++ )
    {
        *intit=(pow(*intit,2));
        sumsq+=*intit;
    }
    sd = sumsq/(deviations.size()-1);
    sd= sqrt(sd);

    statsStruct returnValue;

    returnValue.sd=sd;
    returnValue.sum=sumsq;
    returnValue.mean=mean;
    return returnValue;
}

vector<uRegion> bedgraphToRegions(uTagsExperiment ourExp,const vector<bedScores> vecBedscore, int size, string chromName)
{

    uRegion tempRegion;
    tempRegion.setChr(chromName);
    vector<bedScores>::const_iterator bedit;
    // vector<float>::iterator floatit;
    vector<int>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    vector<uRegion> regionsKept, allregions;

    for ( bedit=vecBedscore.begin(); bedit != vecBedscore.end(); bedit++ )
    {
        {
            tempRegion.setEnd((bedit->position)+30);
            tempRegion.setStart((bedit->position)-30);
            tempRegion.setScore(bedit->score);
            tempRegion.setCount(RegionTags::getTagCount(tempRegion,&ourExp,OverlapType::OVERLAP_CENTRE));
            allregions.push_back(tempRegion);
        }
    }
    return allregions;
}

//Give a bin with list of scores and we generate a region for each score. Optional Treshold
vector<uRegion> makeRegionsfromBin(uNucleoBin& ourBin, uTagsExperiment& ourExp, float threshold)
{

    typedef vector <float> float_v;
    typedef vector <int> int_v;
    uRegion tempRegion;
    tempRegion.setChr(ourBin.getChromName());
    vector<bedScores>::iterator it;
    vector<float>::iterator floatit;
    vector<int>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    vector<uRegion> allregions;
    cerr << "Entered MakeRegionsfromBin" <<endl;

    try
    {
        float_v* scoreVec= ourBin.getpMapScore();
        int_v* countVec=ourBin.getpMapCount();

        //ourBin.writebedGraph(cout);

        //Iterate through our scores given by our bin format and get the regions densities.
        //Testing a refactor on this using the previously measured bin density
        //Create a list of regions. For each position.

        int l=0;
       // cerr << "Entering for"<<endl;
        for ( floatit=scoreVec->begin(); floatit < scoreVec->end(); floatit++ )
        {
            // cerr << "before if" << l <<endl;
            if (*floatit>threshold)
            {
              //  cerr << "in if with end, start and score as "<<(floatit-scoreVec->begin())+30<< " " <<(floatit-scoreVec->begin())-30<<" " << *floatit  <<endl;

                int end = ((floatit-scoreVec->begin())+REGSIZE);
                int start =((floatit-scoreVec->begin())-REGSIZE);

                if (end>=(int)countVec->size())
                    end=countVec->size();
                if (start < 0)
                    start =0;


                tempRegion.setEnd(end);
                tempRegion.setStart(start);
                tempRegion.setScore(*floatit);
                for (int l=tempRegion.getStart(); l<=tempRegion.getEnd(); l++)
                    tempRegion.setCount(tempRegion.getCount()+countVec->at(l));
                allregions.push_back(tempRegion);
            }
            l++;
        }
    }
    catch(...)
    {
        cerr << "Caught error in makeRegionsfromBin() " << endl;
        throw;
    }
    return allregions;
}




//Normalize the given regions according to there densities and variance, using gaussian Similarity. So at most we double high
void Normalize(float sd, float mean,vector<uRegion> & ourRegions)
{
    vector<uRegion>::iterator regit,ic;

    //We have our mean and density deviation
    float mult;
    float score=0;

    //Ponder every score based on it's density.
    for ( regit=ourRegions.begin(); regit < ourRegions.end(); regit++ )
    {
        score=0;
        if (regit->getCount()<(mean))
        {
            score=((utility::gaussianSim(regit->getCount(),(mean),(sd/2)))*regit->getScore());
            regit->setScore(score);
        }
        else
        {
            mult= (1-utility::gaussianSim(regit->getCount(),(mean),(sd/2)));
            score=(1+mult)*regit->getScore();
            regit->setScore(score);
        }
    }
}

//Find the maxima of our score region.
vector<bedScores> getMaximaFromBedscores(const vector<bedScores> ourBedScores)
{

    //If tempty, return empty.
    if (ourBedScores.size()==0)
        return ourBedScores;

    vector<bedScores> returnBed;
    vector<bedScores>::const_iterator bedsit;
    bedScores temp;
    string chr;
    int curPosition, curEnd,lastPosition, lastEnd, lastValue;
    float curValue;
    bool rising=false;

    bedsit=ourBedScores.begin();
    chr = bedsit->chr;
    lastPosition=bedsit->position;
    lastEnd=bedsit->end;
    lastValue=bedsit->score;

    for((bedsit=ourBedScores.begin()+1); bedsit!=ourBedScores.end(); bedsit++)
    {
        chr = bedsit->chr;
        curPosition=bedsit->position;
        curEnd=bedsit->end;
        curValue=bedsit->score;
        //Not adjacent, reset
        if (curPosition==(lastPosition+1))
        {
            //Last position was smaller, we are looking for maxima so continue
            if (lastValue<curValue)
            {
                rising=true;
                lastValue=curValue;
                lastEnd= curEnd;
                lastPosition=curPosition;
            }
            else
            {
                //Last was a maximal, output
                if (rising)
                {
                    temp.chr=chr;
                    temp.position=lastPosition;
                    temp.end=lastEnd;
                    temp.score=lastValue;
                    returnBed.push_back(temp);
                    // cout << chr << "\t" << lastPosition << "\t" <<  lastEnd<< "\t" << lastValue << endl;
                }
                lastValue=curValue;
                lastEnd= curEnd;
                lastPosition=curPosition;
                rising=false;
            }
        }
        else
        {
            if (rising)
            {
                temp.chr=chr;
                temp.position=lastPosition;
                temp.end=lastEnd;
                temp.score=lastValue;
                returnBed.push_back(temp);
                lastValue=curValue;
                lastEnd= curEnd;
                lastPosition=curPosition;
                rising=false;
            }
            lastPosition= curPosition;
            lastEnd = curEnd;
            lastValue= 0.0f;
        }
    }


    return returnBed;
}

//Filter regions below a certain score
vector<uRegion> filterRegions(vector<uRegion> ourRegions, float threshold)
{

    vector<uRegion>::iterator regit;
    vector<uRegion> returnVec;
    for(regit=ourRegions.begin(); regit!=ourRegions.end(); regit++)
    {

        if (regit->getScore()>=threshold)
        {
            returnVec.push_back(*regit);
        }
    }

    return returnVec;

}

//TODO replace this..
//Single use only, makes a bedGraph vector froms a region by using the middle.
vector<bedScores> getBedFromEnlargedRegions(vector<uRegion> & ourRegions, int enlargeSize)
{

    vector<uRegion>::iterator regit;
    vector<bedScores> returnBed;
    bedScores tempBed;
    for ( regit=ourRegions.begin(); regit < ourRegions.end(); regit++ )
    {
        tempBed.chr=regit->getChr() ;
        tempBed.position=(regit->getStart()+(enlargeSize-1));
        tempBed.end=(regit->getEnd()-enlargeSize);
        tempBed.score=regit->getScore();
        returnBed.push_back(tempBed);
    }

    return returnBed;
}

*/
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "uTags.h"
#include "uRegion.h"
#include "functions.h"

using namespace std;


#define REGSIZE 30

std::vector<bedScores> loadbedGraph(string bedpath)
{
    ifstream inputStream;

    inputStream.open(bedpath.c_str());
    if (inputStream.bad())
    {
        cerr << "WTF??. Error loading in loadbedGraph()";
        abort();
    }

    string lineString,header;
    stringstream Infostream;
    std::getline(inputStream, lineString);
    Infostream.str(lineString);
    Infostream >>header;

    if (header.find("type=bedGraph")==string::npos)
    {
        cerr << "Error, this is not a valid bedgraph file"<<endl;
        abort();
    }
    bedScores tempscore;
    vector<bedScores> returnVector;
    while(inputStream.eof()!=true)
    {

        getline(inputStream, lineString);
        //Reset error flags our nothing work.
        Infostream.clear();
        Infostream.str(lineString);
        Infostream >>tempscore.chr;
        Infostream >>tempscore.position;
        Infostream >>tempscore.end;
        Infostream >>tempscore.score;
        returnVector.push_back(tempscore);
    }
    return returnVector;
}


std::vector<bedScores> loadbedGraph(std::ifstream& inputStream)
{
   // ifstream inputStream;

   // inputStream.open(bedpath.c_str());
    if (inputStream.bad())
    {
        cerr << "WTF??. Error loading in loadbedGraph()";
        abort();
    }

    string lineString,header;
    stringstream Infostream;
    std::getline(inputStream, lineString);
    Infostream.str(lineString);
    Infostream >>header;

    if (header.find("type=bedGraph")==string::npos)
    {
        cerr << "Error, this is not a valid bedgraph file"<<endl;
        abort();
    }
    bedScores tempscore;
    vector<bedScores> returnVector;
    while(inputStream.eof()!=true)
    {

        getline(inputStream, lineString);
        //Reset error flags our nothing work.
        Infostream.clear();
        Infostream.str(lineString);
        Infostream >>tempscore.chr;
        Infostream >>tempscore.position;
        Infostream >>tempscore.end;
        Infostream >>tempscore.score;
        returnVector.push_back(tempscore);
    }
    return returnVector;
}


//Pass a list of maxima to this function and it will select specific peaks within competing regions
//As we only want to remove one potentional nucleosome per position.
vector<uRegion> filterNucleoMaxima(uTagsExperiment* ourTagExp, string chr, vector<bedScores> &vecScores)
{

    uRegion tempRegion;
    tempRegion.setChr(chr);
    vector<bedScores>::iterator it;
    vector<int>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    vector<uRegion> regionsKept, allregions;
    // load our maxima
    cerr << " We have " << vecScores.size() << endl;
    cerr << "Loading our maxima" << endl;
    for ( it=vecScores.begin(); it < vecScores.end(); it++ )
    {

        tempRegion.setEnd(it->position+30);
        tempRegion.setStart(it->position-30);
        tempRegion.setScore(it->score);
        //TODO fix this
        tempRegion.setCount( ourTagExp->getSubsetCount(tempRegion,OverlapType::OVERLAP_CENTRE));
        allregions.push_back(tempRegion);
    }

    cerr << "Checking our regions" << endl;
    for ( regit=allregions.begin(); regit < allregions.end(); regit++ )
    {
        bool ok=true;
        ic= regit;
        //Go backwards until before the region
        while ((ic->getEnd()>=regit->getStart())&&(ic!=allregions.begin()))
        {
            ic--;
        }
        //Go forward until past
        while ((ic->getStart()<=regit->getEnd())&&(ic!=allregions.end()))
        {
            if (ic->getScore()>regit->getScore())
            {
                ok =false;
                break;

            }
            ic++;
        }

        if (ok)
            regionsKept.push_back(*regit);
    }

    return regionsKept;
}


uTagsExperiment decompPass(uTagsChrom* ourTagChrom, vector<uRegion> vecRegions)
{

    uTagsExperiment returnTagsExp;
    uTagsChrom tempChrom;
    uRegion tempRegion;
    vector<uRegion>::iterator regit;
    vector<uRegion> regionsRemoved;
    //For each maxima, create a region, subset the tags, etc.
    for ( regit=vecRegions.begin(); regit < vecRegions.end(); regit++ )
    {
        tempChrom=ourTagChrom->getSubset(regit->getStart(), regit->getEnd(), OverlapType::OVERLAP_CENTRE);

         /**< Increase depth vector for every tag we are using */

        returnTagsExp.combine(tempChrom);


    }


    //Remove all our tags from the original dataSet
    int count=0;
    for ( regit=vecRegions.begin(); regit < vecRegions.end(); regit++ )
    {
        count++;
        ourTagChrom->removeSubset(regit->getStart(), regit->getEnd(),OverlapType::OVERLAP_CENTRE);
        //   cerr << "Completed region " << count << endl;
    }



    return returnTagsExp;
}


void decomp(int argc, char* argv[])
{

    if(argc <=5)
    {
        cerr<<"Program signature is -s <samFile> -c <chromosone> -b <BedGraph Maxima> -o <SamLeft> <SamDecomp>" <<endl;
        return;
    }
    string sampath;
    string chromName;
    string bedpath;//; = argv[4];
    string samLeft;//; = argv[5];
    string samDec;// = argv[6];
    uTagsExperiment ourExp;
    uTagsChrom* pChrom;
    ifstream inputStream;
    ifstream bedStream;
    ofstream ofLeft, ofDec;
    stringstream Infostream;
    vector<bedScores> scoreVec;
    bedScores tempScore;

//Parse parameters
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-s")==0)
        {
            sampath = argv[i + 1];
        }
        else if (strcmp(argv[i],"-c")==0)
        {
            chromName = argv[i + 1];
        }
        else if (strcmp(argv[i],"-b")==0)
        {
            bedpath= argv[i + 1];
        }
        else if(strcmp(argv[i],"-o")==0)
        {
            samLeft = argv[i + 1];
            samDec = argv[i + 2];
        }
    }
//Validate we loaded them all
    if ((sampath.size()==0)||(chromName.size()==0)||(bedpath.size()==0)||(samLeft.size()==0)||(samDec.size()==0))
    {
        cerr<<"Program signature is -s <samFile> -c <chromosone> -b <BedGraph Maxima> -o <SamLeft> <SamDecomp>" <<endl;
        return;
    }

    inputStream.open(sampath.c_str());
    if (inputStream.good())
        ourExp.loadFromSam(inputStream);
    //ourExp.loadFromSam(inputStream);
    else
    {
        cerr << "WTF?? Loading";
        abort();
    }
    cerr << "Loading Maxima file" << endl;
    bedStream.open(bedpath.c_str());
    if (bedStream.good())
    {
        string chr;
        string lineString;
        int start, end;
        float score;
        //Skip header
        bedStream>>chr;
        while(bedStream.eof()!=true)
        {
            getline(bedStream, lineString);
            Infostream.clear();
            Infostream.str(lineString);
            Infostream >> chr;
            Infostream >> start;
            Infostream >> end;
            Infostream >> score;

            tempScore.position=start;
            tempScore.score=score;
            scoreVec.push_back(tempScore);
        }
    }
    else
    {
        cerr << "WTF?? Loading on Bed";
        abort();
    }
    //  uTagsChrom* pChrom;
    // vector<int> densityMap;
    pChrom=ourExp.getpChrom(chromName);
    vector<uRegion> regionVector;
    //createDensityMap(pChrom,densityMap);
    cerr << "Filtering our maxima" << endl;
    regionVector=filterNucleoMaxima(&ourExp, chromName,scoreVec);

    ofLeft.open(samLeft.c_str());
    ofDec.open(samDec.c_str());

    uTagsExperiment expDec;
    cerr << "Decomp our data" << endl;
    expDec=decompPass(pChrom,regionVector);
    expDec.sortData();
    ourExp.writeSamToOutput(ofLeft);
    expDec.writeSamToOutput(ofDec);
}

statsStruct returnSdfromRegions(vector<uRegion> regionVec)
{

    vector<bedScores>::iterator it;
    vector<float>::iterator floatit;
    vector<float>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    float sum, mean, sd,sumsq;
    sum=0;
    sumsq=0;
    //Establish our score passed on density average and sigma



    vector<float> deviations;
    for ( regit=regionVec.begin(); regit < regionVec.end(); regit++ )
    {
        sum +=regit->getCount();
    }

    mean = sum/regionVec.size();

    for ( regit=regionVec.begin(); regit < regionVec.end(); regit++ )
    {
        deviations.push_back(regit->getCount()-mean);
    }
    for (intit=deviations.begin(); intit < deviations.end(); intit++ )
    {
        *intit=(pow(*intit,2));
        sumsq+=*intit;
    }
     cerr << "Region count is " << regionVec.size() <<endl;
      cerr <<"sum  is"<< sum <<endl;
      cerr <<"sumQ  is"<< sumsq <<endl;
    cerr << "deviation size is" << (deviations.size()-1) <<endl;


    sd = sumsq/(deviations.size()-1);
    sd= sqrt(sd);

    statsStruct returnValue;

    returnValue.sd=sd;
    returnValue.sum=sumsq;
    returnValue.mean=mean;
    return returnValue;
}

vector<uRegion> bedgraphToRegions(uTagsExperiment ourExp,const vector<bedScores> vecBedscore, int size, string chromName)
{

    uRegion tempRegion;
    tempRegion.setChr(chromName);
    vector<bedScores>::const_iterator bedit;
    // vector<float>::iterator floatit;
    vector<int>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    vector<uRegion> regionsKept, allregions;

    for ( bedit=vecBedscore.begin(); bedit != vecBedscore.end(); bedit++ )
    {
        {
            tempRegion.setEnd((bedit->position)+30);
            tempRegion.setStart((bedit->position)-30);
            tempRegion.setScore(bedit->score);
            tempRegion.setCount(RegionTags::getTagCount(tempRegion,&ourExp,OverlapType::OVERLAP_CENTRE));
            allregions.push_back(tempRegion);
        }
    }
    return allregions;
}


//Give a bin with list of scores and we generate a region for each score. Optional Treshold

vector<uRegion> makeRegionsfromBin(uNucleoBin& ourBin, uTagsExperiment& ourExp, float threshold)
{

    typedef vector <float> float_v;
    typedef vector <int> int_v;

    vector<bedScores>::iterator it;
    vector<float>::iterator floatit;
    vector<int>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    vector<uRegion> allregions;
    try
    {
        float_v* scoreVec= ourBin.getpMapScore();
        int_v* countVec=ourBin.getpMapCount();

        //ourBin.writebedGraph(cout);

        //Iterate through our scores given by our bin format and get the regions densities.
        //Testing a refactor on this using the previously measured bin density
        //Create a list of regions. For each position.
        //bool checkEmpty=false;
        int l=0;
        // cerr << "Entering for"<<endl;
        for ( floatit=scoreVec->begin(); floatit < scoreVec->end(); floatit++ )
        {
            {
                // cerr << "before if" << l <<endl;
                if (*floatit>threshold)
                {
                    //  cerr << "in if with end, start and score as "<<(floatit-scoreVec->begin())+30<< " " <<(floatit-scoreVec->begin())-30<<" " << *floatit  <<endl;
                    uRegion tempRegion;
                    tempRegion.setChr(ourBin.getChromName());
                    int end = ((floatit-scoreVec->begin())+REGSIZE);
                    int start =((floatit-scoreVec->begin())-REGSIZE);

                    if (end>=(int)countVec->size())
                        end=(int)countVec->size()-1;
                    if (start < 0)
                        start =0;

                    tempRegion.setEnd(end);
                    tempRegion.setStart(start);
                    tempRegion.setScore(*floatit);
                    tempRegion.setCount(0);
                    for (int l=tempRegion.getStart(); l<=tempRegion.getEnd(); l++)
                        tempRegion.setCount(tempRegion.getCount()+countVec->at(l));
                    //  if (countVec->at(tempRegion.getEnd())==0)
                    //          checkEmpty=true;
                    allregions.push_back(tempRegion);
                }
            }
            l++;
        }
    }
    catch(...)
    {
        cerr << "Caught error in makeRegionsfromBin() " << endl;
        throw;
    }
    return allregions;
}

//Normalize the given regions according to there densities and variance, using gaussian Similarity. So at most we double high
void Normalize(float sd, float mean,vector<uRegion> & ourRegions)
{
    vector<uRegion>::iterator regit,ic;

    //We have our mean and density deviation
    float mult;
    float score=0;

    //Ponder every score based on it's density.
    for ( regit=ourRegions.begin(); regit < ourRegions.end(); regit++ )
    {
        score=0;
        if (regit->getCount()<(mean))
        {
            score=((utility::gaussianSim(regit->getCount(),(mean),(sd/2)))*regit->getScore());
            regit->setScore(score);
        }
        else
        {
            mult= (1-utility::gaussianSim(regit->getCount(),(mean),(sd/2)));
            score=(1+mult)*regit->getScore();
            regit->setScore(score);
        }
    }
}

/**< Normalize a given score and return it, based on gaussian similarity,  */
float normScore(float score,float density,float sd, float mean)
{

    float  retscore=((utility::gaussianSim(density,(mean),(sd/2)))*score);

    return retscore;
}





//Find the maxima of our score region.
vector<bedScores> getMaximaFromBedscores(const vector<bedScores> ourBedScores)
{

    //If tempty, return empty.
    if (ourBedScores.size()==0)
        return ourBedScores;

    vector<bedScores> returnBed;
    vector<bedScores>::const_iterator bedsit;
    bedScores temp;
    string chr;
    int curPosition, curEnd,lastPosition, lastEnd, lastValue;
    float curValue;
    bool rising=false;

    bedsit=ourBedScores.begin();
    chr = bedsit->chr;
    lastPosition=bedsit->position;
    lastEnd=bedsit->end;
    lastValue=bedsit->score;

    for((bedsit=ourBedScores.begin()+1); bedsit!=ourBedScores.end(); bedsit++)
    {
        chr = bedsit->chr;
        curPosition=bedsit->position;
        curEnd=bedsit->end;
        curValue=bedsit->score;
        //Not adjacent, reset
        if (curPosition==(lastPosition+1))
        {
            //Last position was smaller, we are looking for maxima so continue
            if (lastValue<curValue)
            {
                rising=true;
                lastValue=curValue;
                lastEnd= curEnd;
                lastPosition=curPosition;
            }
            else
            {
                //Last was a maximal, output
                if (rising)
                {
                    temp.chr=chr;
                    temp.position=lastPosition;
                    temp.end=lastEnd;
                    temp.score=lastValue;
                    returnBed.push_back(temp);
                    // cout << chr << "\t" << lastPosition << "\t" <<  lastEnd<< "\t" << lastValue << endl;
                }
                lastValue=curValue;
                lastEnd= curEnd;
                lastPosition=curPosition;
                rising=false;
            }
        }
        else
        {
            if (rising)
            {
                temp.chr=chr;
                temp.position=lastPosition;
                temp.end=lastEnd;
                temp.score=lastValue;
                returnBed.push_back(temp);
                lastValue=curValue;
                lastEnd= curEnd;
                lastPosition=curPosition;
                rising=false;
            }
            lastPosition= curPosition;
            lastEnd = curEnd;
            lastValue= 0.0f;
        }
    }


    return returnBed;
}

//Filter regions below a certain score
vector<uRegion> filterRegions(vector<uRegion> ourRegions, float threshold)
{

    vector<uRegion>::iterator regit;
    vector<uRegion> returnVec;
    for(regit=ourRegions.begin(); regit!=ourRegions.end(); regit++)
    {

        if (regit->getScore()>=threshold)
        {
            returnVec.push_back(*regit);
        }
    }

    return returnVec;

}

//TODO replace this..
//Single use only, makes a bedGraph vector froms a region by using the middle.
vector<bedScores> getBedFromEnlargedRegions(vector<uRegion> & ourRegions, int enlargeSize)
{

    vector<uRegion>::iterator regit;
    vector<bedScores> returnBed;
    bedScores tempBed;
    for ( regit=ourRegions.begin(); regit < ourRegions.end(); regit++ )
    {
        tempBed.chr=regit->getChr() ;
        tempBed.position=(regit->getStart()+(enlargeSize-1));
        tempBed.end=(regit->getEnd()-enlargeSize);
        tempBed.score=regit->getScore();
        returnBed.push_back(tempBed);
    }

    return returnBed;
}


/**< Return the density of every bin */
uNucleoBin loadDensityFromTags(uTagsExperiment& ourExp, int threshold, string chromName,  bool isComplete)
{

    uNucleoBin ourBin;
    ourBin.setChromName(chromName);
    uTagsChrom* pChrom;

    pChrom=ourExp.getpChrom(ourBin.getChromName());


    ourBin.generateDensityBin(&ourExp,pChrom->getChromSize(),threshold, isComplete);

    return ourBin;
}
