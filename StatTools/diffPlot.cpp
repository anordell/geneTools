#include "diffPlot.h"
#include <tclap/CmdLine.h>
#include "utility.h"
#include <string>
#include <string.h>
#include <fstream>
#include <ostream>
#include "uTags.h"
#include "uRegion.h"
//#include "gnuplot_i.hpp"
using namespace std;

int current(0);
int Increment()
{
    return ++current;
}



void diffPlot(int argc, char* argv[])
{
    try
    {
        TCLAP::CmdLine cmd("Generate distribution of distances", ' ', "0.2");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : diffPlot", false," ","String");
        TCLAP::ValueArg<std::string> regionAarg("i","interval","File containing genomic interval of our region",true,"null","Filepath",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Basic name of output file",true,"null","string",cmd);
        TCLAP::ValueArg<int> extendArg("e","ext","Size we want to extend our intervals to",true,0,"int",cmd);
        TCLAP::ValueArg<int> binArg("b","bin","Size of the bins to compress our signal. must be a divisor of ext",false,0,"int",cmd);

        TCLAP::MultiArg<std::string> nameArg("n","name","Name of each mark associated with a -s parameter.",true,"Mark Name",cmd);
        TCLAP::MultiArg<std::string> samArg("s","sam","Path to an epigenetic Sam file, must be two of them",true,"FilePath",cmd);

//       TCLAP::SwitchArg rawSwitch("r", "raw", "If included, output raw signal per interval",cmd);


        cmd.add(command);
        /**< Parse */
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */
        std::string regionAPath = regionAarg.getValue();

//        bool outputRaw=rawSwitch.getValue();
        auto samPathVector = samArg.getValue();
        auto nameVector = nameArg.getValue();
        auto outputName= outputArg.getValue();
        int extendSize=extendArg.getValue();
        int binSize=binArg.getValue();

        /**< Need same number of ID and files and must be 2 */
        if ((samPathVector.size()!=nameVector.size())||(nameVector.size()!=2))
            throw 15;

        /**< Load two marks. and region*/
        //TODO Refactor this to allow N amount of intervals and marks
        uTagsExperiment markA, markB;
        uRegionExperiment regionA, regionB;
        ifstream samStreamA, samStreamB, regionStreamA,regionStreamB;
        utility::loadStream(samPathVector.at(0), samStreamA);
        utility::loadStream(samPathVector.at(1), samStreamB);
        utility::loadStream(regionAPath, regionStreamA);
        utility::loadStream(regionAPath, regionStreamB);

        markA.loadFromSam(samStreamA);
        markB.loadFromSam(samStreamB);
        regionA.loadFromTabFile(regionStreamA);
        regionB.loadFromTabFile(regionStreamB);
        /**< Finished loading data */

        /**< Set region to given size */
        setRegionsSize(regionA, extendSize);
        setRegionsSize(regionB, extendSize);
        /**< Generate Signal for each Mark */


        regionA.measureDensityOverlap(markA);
        try {
        regionA.generateSignal(markA);
        }
        catch(skipped_elem_throw & e)
        {
            /**< Handle invalid elements given to generateSignal */

           cerr << " WARNING: While generating Signal, various tags or regions have been skipped "<<endl;
           cerr << " Mainly, this happens if tags mapped beyond the reference of if your interval file including regions on chr that are not present in the Sam file"<< endl;

           if (vector<uTags> const * vecU =boost::get_error_info<skipped_tags>(e) )
            {
                if (vecU->size()>0)
                {
                    cerr << "Skipped "<<vecU->size()<< " tags, writing them to file tagErrors.txt" <<endl;
                    ofstream output("tagErrors.txt");
                    for (auto & x: *vecU)
                        x.writeSamToOutput(output);
                }
            }
           if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
           {
                if (vecR->size()>0)
                {
                    cerr << "Skipped "<<vecR->size()<< " regions, writing them to file regionErrors.txt" <<endl;
                    ofstream output("regionErrors.txt");
                    for (auto & x: *vecR)
                        x.writeAll(output);

                    cerr <<"Finished writing" <<endl;
                }
           }

        }

                regionB.measureDensityOverlap(markB);
         try {
                regionB.generateSignal(markB);
          }
           catch(skipped_elem_throw & e)
        {
            /**< Handle invalid elements given to generateSignal */

           cerr << " WARNING: While generating Signal, various tags or regions have been skipped "<<endl;
           cerr << " Mainly, this happens if tags mapped beyond the reference of if your interval file including regions on chr that are not present in the Sam file"<< endl;

           if (vector<uTags> const * vecU =boost::get_error_info<skipped_tags>(e) )
            {
                if (vecU->size()>0)
                {
                    cerr << "Skipped "<<vecU->size()<< " tags, writing them to file tagErrors.txt" <<endl;
                    ofstream output("tagErrors.txt");
                    for (auto & x: *vecU)
                        x.writeSamToOutput(output);
                }
            }
           if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
           {
                if (vecR->size()>0)
                {
                    cerr << "Skipped "<<vecR->size()<< " regions, writing them to file regionErrors.txt" <<endl;
                    ofstream output("regionErrors.txt");
                    for (auto & x: *vecR)
                        x.writeAll(output);
                    cerr <<"Finished writing" <<endl;
                }
           }

        }

        /**< Compare distance between every signal and store data. */
        auto distResult=compareSignals(regionA,regionB);

        ofstream outputStream(outputName);
         ofstream allStream(outputName+"allinfo");
       outputStream << "#Similarity data between"<<nameVector.at(0)<< " and "<< nameVector.at(1) <<endl;
       for(auto elem:distResult){
            outputStream << elem.score<<endl;
            allStream << elem.Elem.getChr()<<"\t"<<elem.Elem.getStart()<<"\t"<<elem.Elem.getEnd()<<"\t"<<elem.score << endl;

        }
    }
    catch(int & e)
    {
        if (e==12)
            utility::stringTocerr("in setRegionSize, extendSize parameter <=0");
    }
}

vector<Elem_score> compareSignals(uRegionExperiment & markA,uRegionExperiment & markB)
{

    vector<Elem_score> returnDistances;
    utility::stringTocerr("Comparing Signals");
    for (auto chrIt =markA.first(); chrIt != markA.last(); chrIt++)
    {
       //  utility::stringTocerr("Getting chrom name :"+chrIt->first);
        auto chrA=chrIt->second;
        auto chrB=markB.getpChrom(chrIt->first);
     //   utility::stringTocerr("A has :"+utility::numberToString(chrA.count()));
      //  utility::stringTocerr("B has :"+utility::numberToString(chrB->count()));
        for (int i=0; i < chrA.count(); i++)
        {
            if ( (chrA.getSite(i).getSignal().size())&&(chrB->getSite(i).getSignal().size() )) {
            /**< Generate our distance matrix */
            Elem_score tempElement;
            tempElement.Elem=chrA.getSite(i);
            tempElement.score=(clustering::align_distance(chrA.getSite(i).getSignal(),chrB->getSite(i).getSignal(),0,chrB->getSite(i).getSignal().size()/10));
            returnDistances.push_back(tempElement);
            }

             /**< Debugging info */
        }
    }
    return returnDistances;
}
/**< Extend or shrink the data to fit a given size */
void setRegionsSize(uRegionExperiment & ourRegionExp, int extendSize)
{

    if (extendSize <=0)
        throw 12;

    ourRegionExp.applyOnSites([&](uRegion & Elem)
    {
        int lenght= Elem.getLenght();
        if (lenght>extendSize)
        {
            int over=Elem.getLenght()-extendSize;
            if (over%2)
                Elem.trimSites(over/2,(over/2+1));
            else
                Elem.trimSites(over/2);
        }
        else if (lenght<extendSize)
        {

            int under=extendSize-Elem.getLenght();
            if (under%2)
                Elem.extendSite(under/2,(under/2+1));
            else
                Elem.extendSite(under/2);

        };

    }  );
}


