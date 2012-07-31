#include <tclap/CmdLine.h>
#include "peakStats.h"
#include "uFormats.h"
#include "uTags.h"
#include "uRegion.h"

#include "utility.h"
#include <time.h>


using namespace std;

void peakStatistics(int argc, char* argv[])
{
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : peakStat", false," ","String");
        TCLAP::ValueArg<std::string> nameArg("f","file","File containing genomic regions",true,"null","string",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Path to output file",false,"","string",cmd);



        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );
        /**< Assign */
        std::string filePath = nameArg.getValue();
        std::string outputPath = outputArg.getValue();



        ifstream fileStream;
        utility::loadStream(filePath, fileStream);

        GenomicFileType fileChoice;
        // Do what you intend.

        fileChoice=GenomicFileType::BED;

        if (outputPath.size()!=0)
        {
            ofstream outputOS(outputPath);

            loadData(fileStream, fileChoice, outputOS);
        }
        else
            loadData(fileStream, fileChoice, cout);



    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (int e)
    {

        if (e==10)
            std::cerr << "Error opening file path"<< std::endl;
    }
}

/**< Load Data and write our peaks straight away. */
void loadData(ifstream& fileStream, const GenomicFileType fileChoice,std::ostream& out)
{

    uGenericExp ourExp;

    if(fileChoice==GenomicFileType::BED)
    {
        ourExp.loadFromTabFile(fileStream);
        //factory::
    }

    vector<float> peakSizes;

    for (auto chromIT= ourExp.first(); chromIT!=ourExp.last(); chromIT++)
    {
        (*chromIT).second.applyOnAllSites( [&] (uGenericNGS Elem)
        {
            peakSizes.push_back(Elem.getLenght());
        } );
    }


    auto qResult = utility::quartilesofVector(peakSizes);

    float mean = utility::getMean(peakSizes);

    auto sd = utility::getSd(peakSizes, mean );


    const string tab="\t";
    /**< Header */
    out << "mean" << tab << "sd" << tab << "q1" << tab <<"q2" <<tab<<"q3" <<endl;
    out << mean << tab << sd << tab << qResult.at(0) << tab <<qResult.at(1)<<tab<<qResult.at(2) <<endl;

}




