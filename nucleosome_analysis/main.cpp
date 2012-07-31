#include "uTags.h"
#include "uRegion.h"
#include "uFunctions.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

void establishRegionFromBedGraph(ifstream &bedGraph, uTagsExperiment* ourExperiment, std::vector<uRegion> &regionVec, int leftShift, int rightShift);
void generateSignals(ofstream &out, float lowThreshold, float highThreshold, int mincount, std::vector<uRegion> &regionVec);
void writeSignalLenghts(uTagsExperiment* ourExperiment,int threshold,  int bin, string strand);
void writeDistoGraph(uTagsExperiment* ourExperiment,int threshold,  int bin, string strand);
bool decode_comline(int argc,char **argv);

//void writeRegionInfo(uTagsExperiment* ourExperiment, std::vector<uRegion> &regionVec, std::string format=="TAB");
int main(int argc, char* argv[])
{
    ifstream bedGraph, samFile;
    ofstream outputA;
    string chrom;
    vector<float> ourScore;
    uTagsExperiment expTag;
    std::vector<uRegion> regionVec;


    if (argc<2)
    {
        cerr <<" No argument, quitting"<<endl;
        return 0;

    }


    string firstArg =argv[1];

    if (firstArg=="WriteRegionSignalsFromBedGraph")
    {

        if (argc <7)
        {

            cerr << "Arguments are <inputSam.Sam> <InputRegion.bedgraph> <LeftShift> <RightShift> <outputFile> ";
            return 0;
        }
        samFile.open(argv[2]);
        bedGraph.open(argv[3]);

        int leftshift, rightshift;
        leftshift= atoi(argv[4]);
        rightshift= atoi(argv[5]);

        outputA.open(argv[6]);


        //atof(argv[3]), atof(argv[4]),atoi(argv[6])
        cerr <<"Trying to load " << argv[2] << endl;
        expTag.loadFromSam(samFile, true);
        establishRegionFromBedGraph(bedGraph, &expTag,regionVec, leftshift, rightshift);

        int count;
        float low, high;
        while(true)
        {
            cout << "Input low, high, count" << endl;
            cin >> low;
            cin >> high;
            cin >>count;
            generateSignals(outputA,low, high, count,regionVec);
        }
    }

    if (firstArg=="AVGRegionSignal")
    {

        if (argc <5)
        {

            cerr << "Arguments are <inputSam.Sam> <InputRegion.bedgraph> <outputFile>";
            return 0;
        }
        samFile.open(argv[2]);
        bedGraph.open(argv[3]);

        outputA.open(argv[4]);
        //atof(argv[3]), atof(argv[4]),atoi(argv[6])
        cerr <<"Trying to load " << argv[2] << endl;
        expTag.loadFromSam(samFile, true);
        establishRegionFromBedGraph(bedGraph, &expTag,regionVec, 300, 300);

        int count;
        float low, high;
        while(true)
        {
            cout << "Input low, high, count" << endl;
            cin >> low;
            cin >> high;
            cin >>count;
            generateSignals(outputA,low, high, count,regionVec);
        }
    }


    if (firstArg=="FragmentSize")
    {

        if (argc <6)
        {
            cerr << "Arguments are <inputSam.Sam> <Bin Size> <strand [+,-,x]> <threshold>";
            return 0;
        }

        samFile.open(argv[2]);
        expTag.loadFromSam(samFile, true);

        std::string strand = argv[4];

        if ((strand!="+")&&(strand!="-")&&(strand!="x"))
        {
            cerr << "strand parameter must be +, - or x" <<endl;
            abort();
        }
        cerr <<"Trying to measure Fragment duplication" <<endl;
        cerr << "threshold is" << atoi(argv[5])<<endl;

        writeSignalLenghts(&expTag,atoi(argv[5]), atoi(argv[3]), strand);
    }
    if (firstArg=="RegionInfo")
    {
        if (argc <6)
        {
            cerr << "Arguments are <inputSam.Sam> <InputRegion.bedgraph> <Left Extend> <Right Extend> [optionalOutputFile]";
            return 0;
        }

        samFile.open(argv[2]);
        bedGraph.open(argv[3]);
        int leftshift, rightshift;
        leftshift= atoi(argv[4]);
        rightshift= atoi(argv[5]);

        if (argc>=7)
            outputA.open(argv[6]);
        //atof(argv[3]), atof(argv[4]),atoi(argv[6])
        cerr <<"Trying to load " << argv[2] << endl;
        expTag.loadFromSam(samFile, true);
        establishRegionFromBedGraph(bedGraph, &expTag,regionVec, leftshift, rightshift);
        cerr <<"Trying to establish our regions" << endl;
        for (unsigned int k=0; k < regionVec.size(); k++)
        {
            if (outputA.is_open())
                regionVec.at(k).writeRegion(outputA);
            else
                regionVec.at(k).writeRegion(cout);
        }
    }

    if (firstArg=="Phasogram")
    {

        if (argc <6)
        {
            cerr << "Arguments are <inputSam.Sam> <Bin Size> <strand [+,-,x]> <threshold>";
            return 0;
        }

        samFile.open(argv[2]);
        expTag.loadFromSam(samFile, true);

        std::string strand = argv[4];

        if ((strand!="+")&&(strand!="-")&&(strand!="x"))
        {
            cerr << "strand parameter must be +, - or x" <<endl;
            abort();
        }
        cerr <<"Trying to measure Fragment duplication" <<endl;
        cerr << "threshold is" << atoi(argv[5])<<endl;

        writeDistoGraph(&expTag,atoi(argv[5]), atoi(argv[3]), strand);
    }



}

//Currently non functional.
//
void writeDistoGraph(uTagsExperiment* ourExperiment,int threshold,  int bin, string strand)
{
    uTagsChrom* ptempChrom;
    ptempChrom=ourExperiment->getpChrom("chr21");
//Sort our data, just in case
    //cerr << "Sorting our data. might take a while" << endl;
    // ourExperiment->sortData();
//Every tag
    uTags tempTag;
    vector<int> result;
    vector<int> notresult;
    int prePos=0;
    vector<uTags> ourTags;
    vector<uTags> otherSignal;

    for (int i=0; i< ptempChrom->count(); i++)
    {
        if (i%100000==0)
            cerr<<"Working on tag " <<i << " of " << ptempChrom->count()<<endl;
        if (ptempChrom->getTag(i).getStrand()=='+')
        {
            tempTag=ptempChrom->getTag(i);
            //If same position, keep adding
            if ((tempTag.getStart()-prePos)<=10)
            {
                prePos=tempTag.getStart();
                ourTags.push_back(tempTag);
            }//If not, validate and restart
            else
            {
                //If above Treshold, the last position of our
                if (ourTags.size()>=threshold)
                {
                    int dist=0;
                    //for (int k=0; k<ourTags.size(); k++){
                    dist= tempTag.getStart()-prePos;
                    if (dist+1>result.size())
                        result.resize(dist+1);
                    //Increment our counter
                    result.at(dist)++;
                    //  }
                }/*
                else
                {
                     for (int k=0; k<ourTags.size(); k++){
                        if (ourTags.at(k).getLenght()+1>notresult.size())
                            notresult.resize(ourTags.at(k).getLenght()+1);
                        //Increment our counter
                        notresult.at(ourTags.at(k).getLenght())++;
                    }


                }*/
                //Clean our data
                prePos= tempTag.getStart();
                ourTags.clear();
                ourTags.push_back(tempTag);
            }

        }
    }

    cout <<"Signal for tags with at least " <<  threshold << " in the same region" << endl;

    if (result.size()>2000)
        result.resize(2000);
    for (int i=0; i< result.size(); i++)
        cout << i << "\t";

    cout <<endl;
    for (int i=0; i< result.size(); i++)
        cout << result.at(i) << "\t";

    cout <<endl;
    /*
     cout <<"Signal of the rest of the tags" <<endl;
        for (int i=0; i< notresult.size(); i++)
         cout << i << "\t";

     cout <<endl;
     for (int i=0; i< notresult.size(); i++)
         cout << notresult.at(i) << "\t";

    cout <<endl;*/
}



void writeSignalLenghts(uTagsExperiment* ourExperiment,int threshold,  int bin, string strand)
{

    uTagsChrom* ptempChrom;
    ptempChrom=ourExperiment->getpChrom("chr21");
//Sort our data, just in case
    cerr << "Sorting our data. might take a while" << endl;
    // ourExperiment->sortData();
//Every tag
    uTags tempTag;
    vector<int> result;
    vector<int> notresult;
    int prePos=0;
    vector<uTags> ourTags;
    vector<uTags> otherSignal;

    for (int i=0; i< ptempChrom->count(); i++)
    {
        if (i%100000==0)
            cerr<<"Working on tag " <<i << " of " << ptempChrom->count()<<endl;
        if (ptempChrom->getTag(i).getStrand()=='+')
        {
            tempTag=ptempChrom->getTag(i);

            //If same position, keep adding
            if (tempTag.getStart()==prePos)
            {
                ourTags.push_back(tempTag);
            }//If not, validate and restart
            else
            {
                if (ourTags.size()>=threshold)
                {
                    for (int k=0; k<ourTags.size(); k++)
                    {
                        if (ourTags.at(k).getLenght()+1>result.size())
                            result.resize(ourTags.at(k).getLenght()+1);
                        //Increment our counter
                        result.at(ourTags.at(k).getLenght())++;
                    }
                }
                else
                {
                    for (int k=0; k<ourTags.size(); k++)
                    {
                        if (ourTags.at(k).getLenght()+1>notresult.size())
                            notresult.resize(ourTags.at(k).getLenght()+1);
                        //Increment our counter
                        notresult.at(ourTags.at(k).getLenght())++;
                    }


                }
                //Clean our data
                prePos= tempTag.getStart();
                ourTags.clear();
                ourTags.push_back(tempTag);
            }

        }
    }

    cout <<"Signal for tags with at least " <<  threshold << " starting at same position" << endl;
    for (int i=0; i< result.size(); i++)
        cout << i << "\t";

    cout <<endl;
    for (int i=0; i< result.size(); i++)
        cout << result.at(i) << "\t";

    cout <<endl;
    cout <<"Signal of the rest of the tags" <<endl;
    for (int i=0; i< notresult.size(); i++)
        cout << i << "\t";

    cout <<endl;
    for (int i=0; i< notresult.size(); i++)
        cout << notresult.at(i) << "\t";

    cout <<endl;
}

//Load our region file
void establishRegionFromBedGraph(ifstream &bedGraph, uTagsExperiment* ourExperiment, std::vector<uRegion> &regionVec, int leftShift, int rightShift)
{
    int start, stop;
    string chrom;
    float score;
    uRegion tempRegion;

//Cut header
    bedGraph >>chrom;
    cerr <<"Started loading our regions" <<endl;
    while (bedGraph.eof()==false)
    {

        bedGraph >> chrom;
        bedGraph >> start;
        bedGraph >> stop;
        bedGraph >> score;

        tempRegion.setChr(chrom);
        tempRegion.setStart(start-leftShift);
        tempRegion.setEnd(stop+rightShift);
        tempRegion.setScore( score);
        regionVec.push_back(tempRegion);
    }
    int numberRegion=0;

    cerr <<"Finished loading our regions" <<endl;
    for (unsigned i=0; i< regionVec.size(); i++)
    {


        RegionTags::getSignal(&regionVec.at(i), ourExperiment, true);
        string chromtest = regionVec.at(i).getChr();
        regionVec.at(i).setCount( RegionTags::getTagCount(regionVec.at(i).getChr(),(regionVec.at(i).getStart()),(regionVec.at(i).getEnd()),ourExperiment));
        numberRegion++;

        if(numberRegion%20000==0)
            cerr <<"Done "<< numberRegion << "regions on "<< regionVec.size()<<endl;

    }
    cerr <<" Finished Regions" << endl;

}

void generateSignals(ofstream &out, float lowThreshold, float highThreshold, int mincount, std::vector<uRegion> &regionVec)
{
    int counta=0, countSkip=0;
    vector<float> SignalAVG, tempVec;

    SignalAVG.resize(regionVec.at(0).getLenght());
    for (unsigned int i=0; i< regionVec.size(); i++)
    {
        if ((regionVec.at(i).getScore()>lowThreshold)&&(regionVec.at(i).getScore()<highThreshold))
        {
            if (regionVec.at(i).getCount()>mincount)
            {
                counta++;
                for(int k=0; k<regionVec.at(i).getLenght(); k++)
                {
                    tempVec=regionVec.at(i).getSignal();
                    if (k>= SignalAVG.size())
                        cerr << "Error at position / Item " << i << " " << k << endl;
                    SignalAVG.at(k)+=tempVec.at(k);
                }
            }
            else
                countSkip++;
        }
    }
    out << "Based on " << counta << " regions. Low :  " << lowThreshold << " High : "<< highThreshold<<  endl;
    for(unsigned int k=0; k<SignalAVG.size(); k++)
    {
        SignalAVG.at(k)=((SignalAVG.at(k))/counta);
        out << SignalAVG.at(k) <<"\t";
    }
    out <<endl;


    cerr << "Signal from " << counta << " regions" <<endl;
    cerr << "We skipped " << countSkip << " regions" << endl;
}


