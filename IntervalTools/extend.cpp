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
#include "uInterval.h"

using namespace std;
using namespace TCLAP;

void extend(int argc, char **argv){
 try
    {
        TCLAP::CmdLine cmd("Extend every line of an interval file while keeping extra information", ' ', "0.3");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Must be extend", true," ","string");
        TCLAP::ValueArg<std::string> inputArg("f","file","Path to interval file",true,"","string",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Path to output interval file",false,"","string",cmd);
        TCLAP::ValueArg<int> extendArg("s","size","Extend size",true,0,"int",cmd);


        TCLAP::SwitchArg isLeftExtend("l","left","Extend towards - (leftwise)", false);
        TCLAP::SwitchArg isRightExtend("r","right","Extend towards + (righttwise)", false);
        TCLAP::SwitchArg isBothExtend("b","bi","Extend bi-directional", false);
        vector<Arg*> xorList = {&isLeftExtend,&isRightExtend,&isBothExtend};

        /**< Left or right */
        cmd.xorAdd( xorList );
        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */

        auto filePath = inputArg.getValue();
        auto interValData= loadData(filePath);
        long int extendSize = extendArg.getValue();
        interValData= extendData(interValData,ExtendType::BI, extendSize);

        interValData.writeAsTabFile(cout);
        /**< Need same number of ID and files */
       // if (fileVector.size()!=IDVector.size())
      //      throw 10;
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

uIntervalExperiment loadData(std::string path){
    try{
        uIntervalExperiment returnExp;
        ifstream dataStream;
        utility::loadStream(path, dataStream);
        returnExp.loadFromTabFile(dataStream);
        return returnExp;
    }
    catch(...){
    throw;
    }
}

uIntervalExperiment extendData(uIntervalExperiment intervalData, ExtendType parType, unsigned long int pShiftSize)
{
    intervalData.applyOnSites([&](uInterval& element){
                               switch(parType){
                                case  ExtendType::BI:
                                        element.extendSite(pShiftSize);
                                        break;
                                case  ExtendType::LEFT:
                                        element.extendSite(pShiftSize,0);
                                    break;
                                case  ExtendType::RIGHT:
                                        element.extendSite(0,pShiftSize);
                                    break;
                              }
    }
);
return intervalData;
}

