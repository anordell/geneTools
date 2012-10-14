#include <tclap/CmdLine.h>
#include "histogramFromSam.h"
#include "utility/utility.h"
#include <string>
#include <string.h>
#include <fstream>
#include <ostream>
#include "uTags.h"
#include "omp.h"


#define BIGNUMBER 999999
#define SMALLNUMBER -999999
using namespace std;
using namespace NGS;
/**< Load a sam File */
void histogramFromSam(int argc, char* argv[])
{
    try
    {
        TCLAP::CmdLine cmd("Output data needed to make a histogram on bin density", ' ', "0.3");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : histogramFromSam", false," ","String");
        TCLAP::ValueArg<std::string> outputArg("o","output","Basic name of output file",true,"null","output filename",cmd);
        TCLAP::ValueArg<int> binArg("b","b","Size of the bins, must be above 0",true,0,"bin size",cmd);
        TCLAP::ValueArg<std::string> samArg("s","sam","Path to a Sam file",true,"null", "FilePath",cmd);
        TCLAP::ValueArg<int> upperArg("u","upp","Upper bound to keep data",false,0,"left bound",cmd);
        TCLAP::ValueArg<int> lowerArg("l","low","Lower bound to keep data",false,0,"right bound",cmd);


        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */

        string outputpath = outputArg.getValue();
        string samPath = samArg.getValue();
        int binSize = binArg.getValue();
        int upperBound,lowerBound;

        if (upperArg.getValue()!=0)
            upperBound=upperArg.getValue();
        else
            upperBound=BIGNUMBER;

        if (lowerArg.getValue()!=0)
            lowerBound=lowerArg.getValue();
        else
            lowerBound=SMALLNUMBER;


        cerr << "Bound are "<< lowerBound<<" "<<upperBound <<endl;
        ifstream samStream;
        uTagsExperiment ourExp;

        try {
            utility::loadStream(samPath,samStream);
            ourExp.loadFromSam(samStream);
        }
        catch(...)
        {
            cerr << "Failed while loading" <<endl;
            throw;
        }
       map<int, int> countMap;
       map<string, map<int,int>> chrMap;
       cerr << "Sorting" <<endl;
       ourExp.sortData();
       cerr <<"EndSort"<<endl;
       cerr <<"Subset"<<endl;


        cerr << "Starting pragma" <<endl;

        for(auto itExp = ourExp.begin(); itExp!=ourExp.end(); itExp++)
        {
//                int curPos=0;
//                int pass=0;


            int chromLimit= itExp->second.getChromSize()/binSize;
            chromLimit--;
            // #pragma omp parallel for
            for(int i=0; i<chromLimit; i++)
            {
                int lower =i*binSize;
                int upper =(((i+1)*binSize)-1);
                int curCount = itExp->second.getSubsetCount(lower,upper);
                chrMap[itExp->second.getChr()][curCount]++;
                countMap[curCount]++;
            }

               /* while (curPos < itExp->second.getChromSize() ){
                    int lower = (curPos+1);
                    int upper=(curPos+binSize);
                    int curCount = itExp->second.getSubsetCount(lower,upper);
                    chrMap[itExp->second.getChr()][curCount]++;
                    countMap[curCount]++;
                    curPos=upper;
                }*/
        }
        cerr <<"End Subset"<<endl;
        ofstream outputStream(outputpath+".binCount.log");
        ofstream outputStreamOutlier(outputpath+".outliers.log");
        ofstream outputStreamQuarts(outputpath+".quartiles.log");
        ofstream outputStreamPersonal(outputpath+".personal.log");
       // ofstream outputStreamQuarts(outputpath+"quartiles.log");

        for(auto x:countMap){
            outputStream << x.first << " " << x.second << endl;
        }

        vector<float> valuevector,partialVector;
        auto  it = valuevector.begin();

        for(auto x:chrMap){
                valuevector.clear();
                partialVector.clear();
                outputStream << x.first << endl;
                for(auto y:x.second){
                    if (y.first>0){
                        valuevector.insert(it,y.first,y.second);
                        it = valuevector.begin();
                        outputStream << y.first << " " << y.second << endl;
                    }

                }

                cerr << "Internal" <<endl;
                sort(valuevector.begin(),valuevector.end());
                auto quarts= utility::quartilesofVector(valuevector);
                auto q1=quarts.at(0);
                auto med=quarts.at(1);
                auto q3=quarts.at(2);
                cerr << "Internal 1" <<endl;
                /**< Inter Quartile range */
                auto IQR=q3-q1;
                auto LF = q1 -(1.5*IQR);
                auto UF=  q3 +(1.5*IQR);
                cerr << "Averaging" <<endl<<endl;


                float sum=0, avg=0;
                if (valuevector.size()){
                    sum=std::accumulate(valuevector.begin(),valuevector.end(),0);
                    avg = sum/valuevector.size();
                }
             //   auto partialIt=;


                auto pred = [&] (const float &value) { return ((value >=LF)&&( value<=UF)) ;};
                cerr << "For chrom"<<x.first <<endl;
                cerr << "Frontiers are LF: " << LF <<" UF: "<<UF<<endl;
                 cerr << "There are " << valuevector.size() << " to copy from " <<endl;
                std::copy_if(valuevector.begin(),valuevector.end(), std::back_inserter(partialVector), pred);
                cerr << "Copied " << partialVector.size() << " from " <<valuevector.size()<<endl;
                cerr << "Averaging on non-outliers with " <<partialVector.size() <<" elements."<<endl;

                float sumPartial=0, avgPartial=0;
                if (partialVector.size()){
                    sumPartial=std::accumulate(partialVector.begin(),partialVector.end(),0);
                    avgPartial = sumPartial/partialVector.size();
                }
                outputStreamQuarts << "For chrom"<<x.first <<endl;
                outputStreamQuarts << "Q1 " << q1 <<" ";
                outputStreamQuarts << "med " << med <<" ";
                outputStreamQuarts << "Q3 " << q3 <<" ";
                outputStreamQuarts << "LF " << LF <<" ";
                outputStreamQuarts << "UF " << UF <<endl;
                outputStreamQuarts << "AVG: " << avg <<" ";
                outputStreamQuarts << "AVG no outliers " << avgPartial <<endl;
                outputStreamQuarts << "SD: " << utility::getSd(valuevector,avg) <<" ";
                outputStreamQuarts << "SD no outliers: " << utility::getSd(partialVector,avgPartial) <<endl;

            //    for(auto y:x.second){
            //        if ((y.first <LF)||(y.first >UF))
            //            outputStreamQuarts << "Outlier value and count is" << y.first << " " << y.second <<endl;
            //    }
                int curPos=0;
                cerr <<"Writing outlier and personal limits" <<endl;
                while (curPos < ourExp.getChrSize(x.first) ){
                    int lower = (curPos+1);
                    int upper=(curPos+binSize);
                    int curCount = ourExp.getSubsetCount(x.first,lower,upper);
                    if (curCount>UF)
                        outputStreamOutlier<<x.first<<"\t"<<lower<<"\t"<<upper << "\t" <<curCount<<endl;

                    if ((upperBound<BIGNUMBER)||(lowerBound>SMALLNUMBER))
                    {
                          if ((curCount<lowerBound)||(curCount>upperBound)){
                              outputStreamPersonal<<x.first<<"\t"<<lower<<"\t"<<upper << "\t" <<curCount<<endl;
                        }
                    }
                    countMap[curCount]++;
                    curPos=upper;
                }
            }

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (ugene_exception_base &e)
    {
        if (std::string const * ste =boost::get_error_info<string_error>(e) )
        {
            cerr << "Trace of crash:"<<endl;
            cerr << *ste;
        }

    }
    catch(exception &e){

        utility::stringTocerr("An error we do not handle???, outputing what()");
        utility::stringTocerr(e.what());
    }
    catch(...){

        utility::stringTocerr("An error we do not handle???");
       // utility::stringTocerr(e.what());
    }
}






