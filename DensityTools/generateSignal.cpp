#include <tclap/CmdLine.h>
#include "generateSignal.h"
#include "utility.h"
#include <string>
#include <string.h>
#include <fstream>
#include <ostream>
#include "uTags.h"
#include "uRegion.h"
#include "gnuplot_i.hpp"

//class Gnuplot;
using namespace std;

int current(0);
int Increment()
{
    return ++current;
}

void generateSignal(int argc, char* argv[])
{
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.4.4");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : peakStat", false," ","String");
        TCLAP::ValueArg<std::string> regionAarg("i","interval","File containing genomic interval regions",true,"null","Filepath",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Basic name of output file",true,"null","string",cmd);
        TCLAP::ValueArg<int> extendArg("e","ext","Size we want to extend our intervals to",true,0,"int",cmd);
        TCLAP::ValueArg<int> binArg("b","b","Size of the bins to compress our signal. must be a divisor of ext",false,0,"int",cmd);

        TCLAP::MultiArg<std::string> nameArg("n","name","Name of each mark associated with a -s parameter.",true,"Mark Name",cmd);
        TCLAP::MultiArg<std::string> samArg("f","sam","Path to an epigenetic Sam or Wig(density) file",true,"FilePath",cmd);

        TCLAP::SwitchArg rawSwitch("r", "raw", "If included, output raw signal per interval",cmd);
        TCLAP::SwitchArg wiggleswitch("w", "wig", "Temporary, include if input files are in wiggle format",cmd);

        cmd.add(command);
        /**< Parse */
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */

        std::string regionAPath = regionAarg.getValue();

        bool outputRaw=rawSwitch.getValue();
        bool wigSwitch= wiggleswitch.getValue();
        auto samPathVector = samArg.getValue();
        auto nameVector = nameArg.getValue();
        auto outputName= outputArg.getValue();
        /**< Need same number of ID and files */
        if (samPathVector.size()!=nameVector.size())
            throw 15;

        int extendSize=extendArg.getValue();
        int binSize=binArg.getValue();
        vector<vector<float>> fullSignal;
        vector<vector<float>> normSignal;
        vector<vector<float>> sdSignal;
        /**< Used only to provided a "position" to each data point */
        vector<float> X_Vector;
        X_Vector.resize(extendSize);
        generate (X_Vector.begin(), X_Vector.end(), Increment);

        /**< For each Mark, do Signal. */
        for (int i=0; i< (int)samPathVector.size(); i++ )
        {

            {
                ifstream samStream;
                //  utility::stringTocerr("Before Load Sam");
                utility::stringTocerr("Opening Sam Stream");
                utility::loadStream(samPathVector.at(i), samStream);
                string markName= nameVector.at(i);
                ifstream regionAStream;
                utility::stringTocerr("Opening interval Stream");
                utility::loadStream(regionAPath, regionAStream);

                /**< Overlap and get density */
                uTagsExperiment tagTreat;
                uRegionExperiment wigTreat;
                utility::stringTocerr("Loading Sam/Wig");

                if (!(wigSwitch)){
                    tagTreat.loadFromSam(samStream);
                     tagTreat.sortData();}
                else
                {
                    wigTreat.loadFromWig(samStream);
                    wigTreat.sortData();}

                uRegionExperiment regionExpCtrl;
                regionExpCtrl.loadFromTabFile(regionAStream);
                utility::stringTocerr("Loaded Interval");

                /**< Normalise our sites */
                setRegionsSize(regionExpCtrl, extendSize);

                utility::stringTocerr("Measuring Density");

                regionExpCtrl.measureDensityOverlap(tagTreat);
                utility::stringTocerr("Generating Signal");

                try
                {
                if (!(wigSwitch)){
                    regionExpCtrl.generateSignal(tagTreat);
                    }
                else
                    {
                    regionExpCtrl.generateSignal(wigTreat);
                   }
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
                            ofstream output(outputName+markName+"tagErrors.txt");
                            for (auto & x: *vecU)
                                x.writeSamToOutput(output);
                        }
                    }
                    if (vector<uRegion> const * vecR =boost::get_error_info<skipped_regions>(e) )
                    {
                        if (vecR->size()>0)
                        {
                            cerr << "Skipped "<<vecR->size()<< " regions, writing them to file regionErrors.txt" <<endl;
                            ofstream output(outputName+markName+"regionErrors.txt");
                            for (auto & x: *vecR)
                                x.writeAll(output);

                            cerr <<"Finished writing" <<endl;
                        }
                    }

                }
                utility::stringTocerr("Averaging Signal");
                auto tempAVG=move(getAvgSignal(regionExpCtrl,extendSize));
                utility::stringTocerr("Normalizing");

                   vector<float> tempRPM;
               if (!(wigSwitch))
                    tempRPM= move(normRPM(tempAVG,tagTreat));
                /**< Bin our data is needed */

                if (binSize)
                {
                    utility::stringTocerr("Binning");
                    tempAVG=binSignal(tempAVG,binSize);
                    if (!(wigSwitch))
                       tempRPM=binSignal(tempRPM,binSize);
                }
                utility::stringTocerr("Saving Full");
                fullSignal.push_back(tempAVG);
               if (!(wigSwitch)){
                utility::stringTocerr("Saving RPM");
                normSignal.push_back(tempRPM); }
                // sdSignal.push_back(getSDSignal(regionExpCtrl,extendSize));

                /**< Output raw data if switch is on */
                if (outputRaw)
                {
                    cerr << "Output raw" <<endl;
                    ofstream rawOutput(outputName+nameVector.at(i)+".raw.txt");
                    regionExpCtrl.writeSignal(rawOutput);
                }

            }
        }
        ofstream SUMOUT(outputName+".sum.txt");
        ofstream RPM(outputName+".rpm.txt");
        ofstream SUMOUTOrtho(outputName+".orthosum.txt");
        ofstream RPMOrtho(outputName+".orthorpm.txt");

        /**< GnuPlot output */
        utility::stringTocerr("Writing first Ps");
        Gnuplot g1("lines");
        g1.set_style("points");
        g1.savetops(outputName+".sum");
        /**< Output to PS full signal in one graph */
        for (int i=0; i<(int)fullSignal.size() ; i++)
        {
            g1.plot_x(fullSignal.at(i),nameVector.at(i));
        }
        /**< Output to PS full signal individually */
        for (int i=0; i<(int)fullSignal.size() ; i++)
        {
            g1.reset_plot();
            // g1.set_style("lines");
            g1.plot_x(fullSignal.at(i),nameVector.at(i));
        }
        utility::stringTocerr("Writing second Ps");
        Gnuplot g2("lines");
        g2.set_style("points");
        g2.savetops(outputName+".rpm");
        /**< Output to PS full signal in one graph */
        for (int i=0; i<(int)normSignal.size() ; i++)
        {
            g2.plot_x(normSignal.at(i),nameVector.at(i));
        }
        /**< Output to PS full signal individually */
        for (int i=0; i<(int)normSignal.size() ; i++)
        {
            g2.reset_plot();
            // g1.set_style("lines");
            g2.plot_x(normSignal.at(i),nameVector.at(i));
        }

        /**< Row */

        for (int i=0; i<(int)fullSignal.size() ; i++)
        {
            /**< Write Summed data file */
            SUMOUT<< nameVector.at(i) <<endl;
            for (auto x: fullSignal.at(i))
                SUMOUT << x << " ";
            SUMOUT <<endl;
            /**< WRite normalized by RPM file */
            RPM<< nameVector.at(i) <<endl;
            for (auto x: normSignal.at(i))
                RPM << x << " ";
            RPM <<endl;
        }
        /**< Column */
        RPMOrtho<<"#";
        SUMOUTOrtho<<"#";
        /**< Write names */
        for (auto name: nameVector)
        {
            RPMOrtho<<name << "\t";
            SUMOUTOrtho<<name << "\t";
        }
        RPMOrtho<<endl;
        SUMOUTOrtho<<endl;
        for (int i=1; i<(int)fullSignal.at(0).size() ; i++)
        {
            for (int j=0; j<(int)fullSignal.size(); j++ )
            {
                SUMOUTOrtho << fullSignal.at(j).at(i) << "\t";
                RPMOrtho << normSignal.at(j).at(i) << "\t";
            }
            RPMOrtho<<endl;
            SUMOUTOrtho<<endl;
        }

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (int & e)
    {
        if (e==15)
            std::cerr << "Please input an equal number of -n and -s parameters"<< std::endl;
        else if (e==20)
            std::cerr << "Failed loading data"<< std::endl;
        else if (e==11)
        {
            std::cerr << "Bin size must be a valid divider for extSize"<< std::endl;
            std::cerr << "Ext ="<< std::endl;
        }

        else
        {
            cerr << "Some kind of crash, ERROR code is " <<e << std::endl;
            throw;
        }

    }
    catch (elem_throw & e)
    {
        const uRegion* errorRegPoint;
        const uTags* errorTagPoint;
        if( ( errorRegPoint =boost::get_error_info<region_error>(e)) )
        {
            cerr <<" We crashed working on this uRegion" <<endl;
            errorRegPoint->debugElem();
        }
        if( ( errorTagPoint=boost::get_error_info<tag_error>(e) ) )
        {
            cerr <<" We crashed working on this uTag" <<endl;
            errorTagPoint->debugElem();
        }

        if (std::string const * ste =boost::get_error_info<string_error>(e) )
            {
            cerr << "Trace of crash" <<endl;
            cerr << *ste;
        }
    }
    catch (std::exception & e)
    {
        utility::stringTocerr("An error we do not handle, outputing what()");
        utility::stringTocerr(e.what());
    }
    catch(...)
    {
         cerr << "should have caught this before? Closing"<<endl;
    }
}


/**< Extend or shrink the data to fit a given size */
void setRegionsSize(uRegionExperiment & ourRegionExp, int extendSize)
{
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


/** \brief Get the Standard deviation of every position of our signal, must be normalized first
 *
 * \param ourRegionExp uRegionExperiment&
 * \param extendSize int
 * \return vector<float>
 *
 */
vector<float> getSDSignal(uRegionExperiment & ourRegionExp, int extendSize)
{
    vector<vector<float>> returnSignal;
    vector<float> returnSdVec;
    returnSignal.resize(extendSize);
    returnSdVec.resize(extendSize);
    for( auto elm: returnSignal)
        elm.resize(ourRegionExp.count());

    ourRegionExp.applyOnSites([&](uRegion & Elem)
    {
        auto tempSignal= Elem.getSignal();

        for (int i=0 ; i< (int)tempSignal.size(); i++)
        {
            returnSignal.at(i).push_back(tempSignal.at(i));
        }
    }  );

    for (int i=0 ; i< (int)returnSignal.size(); i++)
    {
        auto mean =utility::getMean(returnSignal.at(i));
        returnSdVec.at(i)=utility::getSd(returnSignal.at(i),mean);
    }
    return returnSdVec;
}

/** \brief Get Average Signal, so sum/nyumber of regions
 *
 * \param ourRegionExp uRegionExperiment&
 * \param extendSize int
 * \return vector<float>
 *
 */
vector<float> getAvgSignal(uRegionExperiment & ourRegionExp, int extendSize)
{

    vector<float> returnSignal;
    returnSignal.resize(extendSize);
    int emptyCount=0;
    ourRegionExp.applyOnSites([&](uRegion & Elem)
    {

        auto tempSignal= Elem.getSignal();
        if (tempSignal.size()==0)
        {
            emptyCount++;
            emptyCount++;
            cerr<<" Skipping region with no signal. Output: " <<endl;
            Elem.debugElem();
        }
        for (int i=0 ; i< (int)tempSignal.size(); i++)
        {
            returnSignal.at(i)+=tempSignal.at(i);
        }

    }  );
    if (emptyCount)
        cerr <<emptyCount << " regions had no signal. Most likely mapping off reference or on a chr not represented in the Sam files" <<endl;
    for (int i=0 ; i< (int)returnSignal.size(); i++)
    {
        returnSignal.at(i)/=( ourRegionExp.count()-emptyCount);
    }
    cerr<< "We used "<<( ourRegionExp.count()-emptyCount) <<"regions to make this Average Signal. If this number makes no sense, please send a nasty email to Alexei! "  <<endl;

    return returnSignal;
}

/** \brief Normalized a signal by RPM of an experiment
 *
 * \param input_Signal vector<float>
 * \param ourTags const uTagsExperiment&
 * \return vector<float>
 *
 */
vector<float> normRPM(vector<float> input_Signal,const uTagsExperiment & ourTags)
{

    const long int ONE_MILLION=1000000;

    long int tagCount= ourTags.count();

    int normVal= tagCount/ONE_MILLION;

    for (int i=0; i< (int)input_Signal.size(); i++)
    {
        input_Signal.at(i)/=normVal;
    }
    return input_Signal;
}



/** \brief
 *
 * \param input_Signal const vector<float>&
 * \param binSize const int&
 * \return vector<float>
 *
 */
vector<float> binSignal(const vector<float> & input_Signal, const int & binSize)
{

    if (input_Signal.size()%binSize!=0)
        throw 11;

    vector<float> returnSignal;
    returnSignal.resize(input_Signal.size()/binSize);

    /**< For each position in binned vectors, sum the corresponding Bin positions */
    for (int i=0; i< (int)returnSignal.size(); i++)
    {
        for( int x=0; x<binSize; x++ )
            returnSignal.at(i)+=input_Signal.at((i*binSize)+x);
        /**< Divide by Bin Size */
        returnSignal.at(i)/=binSize;
    }

    return returnSignal;
}
