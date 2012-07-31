#include "uInterval.h"
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>
#include "utility.h"
#include <time.h>

using namespace std;

uInterval::uInterval( std::string chr, int start, int end, string extra){

     try {

        setChr(chr);
        setStartEnd(start,end);
        setExtra(extra);
    }
    catch(...){
    cerr <<"Error in constructor of uInterval";
        throw;
    }
 }

uInterval::uInterval(uGenericNGS tag)
{

    setChr(tag.getChr());
    setStartEnd(tag.getStart(),tag.getEnd());
}


void uIntervalExperiment::loadFromTabFile(std::ifstream& stream)
{
    std::string tempString;
    while(!std::getline(stream, tempString).eof())
    {
        addSite(factory::makeIntervalfromTabString(tempString));
    }

}


void uInterval::writeAsTabFile(std::ostream& stream)
{
     stream << chr << "\t" << startPos << "\t" << endPos <<"\t" << m_Extra <<endl;

}


void uIntervalChrom::writeAsTabFile(std::ostream& stream){

 applyOnAllSites(bind2nd(mem_fun_ref(&uInterval::writeAsTabFile), stream));

}


void uIntervalExperiment::writeAsTabFile(std::ostream& stream){

   applyOnAllChroms(bind2nd(mem_fun_ref(&uIntervalChrom::writeAsTabFile), stream));
}

