#include "getDensity.h"

#include <string>
#include <string.h>
#include <fstream>
#include "uTags.h"
#include "utility.h"
#include "uRegion.h"




using namespace std;
/** \brief Taking a tab file and Sam file as input, we return the number of tags from the sam file overlapping our tab regions.
 *
 * \param argc int
 * \param argv[] char*
 * \return void
 *
 */
void getDensity(int argc, char* argv[])
{
    string samPath="";
    string pathname="";
    string outputPath="";
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-f")==0)
        {
            pathname = argv[i + 1];
        }
        else if (strcmp(argv[i],"-s")==0)
        {
            samPath = argv[i + 1];
        }
        else if (strcmp(argv[i],"-o")==0)
        {
            outputPath= argv[i + 1];
        }
    }
    //If we did not input file or chrom name
    if ((pathname.size()==0)||(samPath.size()==0))
    {
        cerr<<"Program signature is density -f <BinPath> -s <SamPath>  -o [OutputPath]";
        return;
    }


    ifstream inputStream(samPath);
    if( !( inputStream))
    {
        cerr <<"Invalid file name for Sam File. File name is "<< pathname <<endl;
        return;
    }
    ifstream tabStream(pathname);
    if( !( tabStream))
    {
        cerr <<"Invalid file name for Tab File. File name is "<< pathname <<endl;
        return;
    }

    /**< Overlap and get density */
    uTagsExperiment tagExp;
    tagExp.loadFromSam(inputStream);

    uRegionExperiment regionExp;
    regionExp.loadFromTabFile(tabStream);

    tagExp.sortData();
    regionExp.measureDensityOverlap(tagExp);
    if (outputPath.size()!=0)
    {
        ofstream outputOS(outputPath);

        regionExp.writeDensityAsTab(outputOS);
    }
    else
        regionExp.writeDensityAsTab(cout);

}
