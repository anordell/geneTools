#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include "utility.h"
#include <time.h>
#include <tclap/CmdLine.h>
#include "extend.h"
#include "uTags.h"
#include "uRegion.h"
#include "subSetSamFromInterval.h"


using namespace std;

void subSetSamFromInterval(int argc, char **argv){
 try
    {
        TCLAP::CmdLine cmd("Extend every line of an interval file while keeping extra information", ' ', "0.3");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<std::string>  command("command", "Must be extend", true," ","string");
        TCLAP::ValueArg<std::string> inputArg("f","file","Path to interval file",true,"","string",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Path to output Sam file",true,"","string",cmd);
        TCLAP::ValueArg<std::string> samArg("s","sam","Path to Sam file",true,"","string",cmd);

        TCLAP::SwitchArg isStrictOverlap("","strict","Only keep tags stricly included in interval", false);


        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */

        auto filePath = inputArg.getValue();
        auto samPath = samArg.getValue();

        auto tagsData=  loadSamData(samPath);
        auto intervalData= loadRegionData(filePath);

       // long int extendSize = extendArg.getValue();
        //interValData= extendData(interValData,ExtendType::BI, extendSize);

//        interValData.writeAsTabFile(cout);
        /**< Need same number of ID and files */
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (...)
    {
        std::cerr << "Some other error" <<endl;
    }
}

uTagsExperiment loadSamData(std::string path){
    try{
        uTagsExperiment returnExp;
        ifstream dataStream;
        utility::loadStream(path, dataStream);
        returnExp.loadFromSam(dataStream);
        return returnExp;
    }
    catch(...){
    throw;
    }
}

uRegionExperiment loadRegionData(std::string path){
    try{
        uRegionExperiment returnExp;
        ifstream dataStream;
        utility::loadStream(path, dataStream);
        returnExp.loadFromTabFile(dataStream);
        return returnExp;
    }
    catch(...){
    throw;
    }
}

void writeSamSubSet(ostream & out,const uRegionExperiment & regionExp, uTagsExperiment tagExp ){






}

