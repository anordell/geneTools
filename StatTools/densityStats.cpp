#include <tclap/CmdLine.h>
#include "densityStats.h"
#include "utility.h"
#include <string>
#include <string.h>
#include <fstream>
#include "uTags.h"
#include "uRegion.h"


using namespace std;

void densityStats(int argc, char* argv[])
{
        TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : peakStat", false," ","String");
        TCLAP::ValueArg<std::string> regionArg("f","file","File containing genomic interval regions",true,"null","string",cmd);
        TCLAP::ValueArg<std::string> samArg("s","sam","Single end sam files path",true,"null","string",cmd);
      //  TCLAP::ValueArg<std::string> outputArg("o","output","Path to output file",false,"","string",cmd);
      //  TCLAP::SwitchArg isBedSwitch("b","bed","If a bed file", cmd, false);
        //     TCLAP::SwitchArg isSamSwitch("s","sam","If a Sam file", cmd, false);
        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */
        std::string regionPath = regionArg.getValue();
        std::string samPath = samArg.getValue();

        ifstream samStream;
        utility::loadStream(samPath, samStream);
        ifstream regionStream;
        utility::loadStream(regionPath, regionStream);

    /**< Overlap and get density */
    uTagsExperiment tagExp;
    tagExp.loadFromSam(samStream);
    tagExp.sortData();

    uRegionExperiment regionExp;
    regionExp.loadFromTabFile(regionStream);

    regionExp.measureDensityOverlap(tagExp);
   // regionExp
    /**< Have to measure density beforeHand */
    auto statsResult = getDensityStats(regionExp);
    auto normstatsResult = getNormDensityStats(regionExp );
    auto coverageStatResult = getCoverageDensityStats(regionExp, tagExp);
    writeData(statsResult,cout);
    cout << "Normalised density per peak" <<endl;
    writeData(normstatsResult,cout);
    cout << "Coverage per bp" <<endl;
    writeData(coverageStatResult,cout);

}

StatsStructure getCoverageDensityStats(uRegionExperiment& regionExp, uTagsExperiment& tagExp){

    vector<float> densityVector;
    int skipped=0;
   for (auto chromIT= regionExp.first(); chromIT!=regionExp.last(); chromIT++)
    {
       {

        vector<long int> densityValues;
        auto pTag= tagExp.getpChrom(chromIT->first);
        densityValues.resize(pTag->getChromSize());
        cerr << "Resized to " <<pTag->getChromSize() << " of chrom " <<chromIT->first <<endl;
        //  auto func = ([&] (uTags Elem){        densityValues.at(Elem.getStart());   } );
       (pTag)->applyOnAllSites([&] (uTags Elem)
        {
            for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
            {
                if (i <(int)densityValues.size())
                    densityValues.at(i)++;
            }
        }

        );
        (*chromIT).second.applyOnAllSites( [&] (uRegion Elem)
        {
            float sum=0;
            float missedBP=0;
            for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
            {
               if (i <(int)densityValues.size())
                sum+= densityValues.at(i);
               else{
                    skipped++;
                    missedBP++;
                    }
            }
            if (missedBP)
              cerr << " Skipping " << missedBP << " Spots with a density of " << sum <<" at pos " <<Elem.getStart() <<" " <<Elem.getEnd() <<endl;
            float effectiveLenght = (float)Elem.getLenght()-missedBP;

            densityVector.push_back(sum/ effectiveLenght) ;
        } );

        }

    }

        cerr << "Skipped " << skipped << endl;

    StatsStructure statResult;

    auto qResult = utility::quartilesofVector(densityVector);

    statResult.q1= qResult.at(0);
    statResult.median =qResult.at(1);
    statResult.q3=qResult.at(2);
    statResult.mean = utility::getMean(densityVector);
    statResult.sd = utility::getSd(densityVector,  statResult.mean );

    return statResult;
}

StatsStructure getNormDensityStats(uRegionExperiment& regionExp){

 try{


    vector<float> densityVector;

    for (auto chromIT= regionExp.first(); chromIT!=regionExp.last(); chromIT++)
    {
        (*chromIT).second.applyOnAllSites( [&] (uRegion Elem)
        {
            densityVector.push_back((float)Elem.getCount()/(float)Elem.getLenght());
        } );
    }

    StatsStructure statResult;

    auto qResult = utility::quartilesofVector(densityVector);

    statResult.q1= qResult.at(0);
    statResult.median =qResult.at(1);
    statResult.q3=qResult.at(2);
    statResult.mean = utility::getMean(densityVector);
    statResult.sd = utility::getSd(densityVector,  statResult.mean );

    return statResult;
    }
    catch(...)
    {
        throw;
    }
}



StatsStructure getDensityStats(uRegionExperiment& regionExp)
{

    try{


    vector<float> densityVector;

    for (auto chromIT= regionExp.first(); chromIT!=regionExp.last(); chromIT++)
    {
        (*chromIT).second.applyOnAllSites( [&] (uRegion Elem)
        {
            densityVector.push_back(Elem.getCount());
        } );
    }


    StatsStructure statResult;

    auto qResult = utility::quartilesofVector(densityVector);

    statResult.q1= qResult.at(0);
    statResult.median =qResult.at(1);
    statResult.q3=qResult.at(2);
    statResult.mean = utility::getMean(densityVector);
    statResult.sd = utility::getSd(densityVector,  statResult.mean );

    return statResult;
    }
    catch(...)
    {
        throw;
    }


}

void writeData(StatsStructure statResult, ostream & out)
{

const string tab="\t";

    out << "mean" << tab << "sd" << tab << "q1" << tab <<"median" <<tab<<"q3" <<endl;
    out << statResult.mean << tab << statResult.sd << tab << statResult.q1 << tab <<statResult.median<<tab<<statResult.q3 <<endl;

}

