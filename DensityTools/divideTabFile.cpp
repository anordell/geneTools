#include "divideTabFile.h"
#include "utility/utility.h"
#include "uRegion.h"
#include <string>
#include <string.h>
#include <fstream>
using namespace std;
using namespace NGS;

/** \brief Divide the given tab delimited file (string, int int ) into a bed file with each region being split into N bins according to parameters
 *
 * \param
 * \param
 * \return
 *
 */

void divideTabFile(int argc, char* argv[])
{

    int size=0;
    string completeType="";
    string pathname="";
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-f")==0)
        {
            pathname = argv[i + 1];
        }
        else if (strcmp(argv[i],"-t")==0)
        {
            completeType = argv[i + 1];
        }
        else if (strcmp(argv[i],"-s")==0)
        {
            size = atoi(argv[i + 1]);
        }
    }
    //If we did not input file or chrom name
    if ((pathname.size()==0)||(size==0))
    {
        cerr<<"Program signature is -f <filepath> -s <Bin Size> -t [default=STRICT; IGNORE/STRICT/EXTEND/ADD]";
        return;
    }

    SplitType type=SplitType::STRICT;
    if (completeType.size()!=0)
    {
        if (completeType=="STRICT")
            type=SplitType::STRICT;
        else if (completeType=="IGNORE")
            type=SplitType::IGNORE;
        else if (completeType=="ADD")
            type=SplitType::ADD;
        else if (completeType=="EXTEND")
            type=SplitType::EXTEND;
        else
        {
            cerr << "Invalid -t parameter. Parameter is "<< completeType <<endl;
            return;
        }
    }


    ifstream inputStream(pathname);
    if( !( inputStream))
    {
        cerr <<"Invalid file name. File name is "<< pathname <<endl;
        return;
    }
    uParser bed3Parser("BED",pathname);
    uWriter bedWriter(&cout,"BED3");
    uRegionExperiment ourExp;
    try
    {
        ourExp.loadWithParser(bed3Parser);
        //  auto ourChrom= ourExp.getpChrom("chr1");

        //  ourChrom->divideItemsIntoBinofSize(size,type);
        //  ourChrom->outputBedFormat(cout);
        ourExp.divideItemsIntoBinofSize(size,type);
        ourExp.writeWithWriter(bedWriter);
    }
    catch(...)
    {
        cerr << "failed. Are your regions smaller then your bin size, or your -t at STRICT ?" <<endl;
    }
}
