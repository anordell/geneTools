#include "densityFromFile.h"
#include "utility.h"
#include <string>
#include <string.h>
#include <fstream>
#include "uTags.h"
using namespace std;
/** \brief Taking a data file as input, this generates a bedgraph style file with the data density per region.
 *
 * \param argc int
 * \param argv[] char*
 * \return void
 *
 */
void densityFromFile(int argc, char* argv[])
{

try{
    string samPath="";
    string pathname="";
    string outputPath="";
    bool bedgraph=false;
    int binSize =0;
    for (int i = 1; i < argc; i++)
    {

        if (strcmp(argv[i],"-s")==0)
        {
            samPath = argv[i + 1];
        }
        else if (strcmp(argv[i],"-o")==0)
        {
            outputPath= argv[i + 1];
        }
        else if (strcmp(argv[i],"-b")==0)
            bedgraph=true;
        else if (strcmp(argv[i],"-bin")==0)
        {
            binSize= atoi(argv[i + 1]);
        }
    }
    //If we did not input file or chrom name
    if ((samPath.size()==0))
    {
        cerr<<"Program signature is densityFromSam -s <SamFile> -o [OutputPath]";
        return;
    }


    ifstream inputStream(samPath);
    if( !( inputStream))
    {
        cerr <<"Invalid file name for Sam File. File name is "<< pathname <<endl;
        return;
    }

    /**< Overlap and get density */
    uTagsExperiment tagExp;
    tagExp.loadFromSam(inputStream);
    cerr <<"Sorting" <<endl;
    tagExp.sortData();


    std::ofstream outputOS;
    if (outputPath.size()!=0)
    {
        outputOS.open(outputPath);
    }

    std::ostream & outFile = ((outputPath.size()!=0) ? outputOS : std::cout);

    cerr << "Starting bin factor" <<endl;
    if (binSize)
    {
        cerr << " Binning" <<endl;
        for (auto chromIT= tagExp.first(); chromIT!=tagExp.last(); chromIT++)
        {
            //auto tagChrom= tagExp.getpChrom("chr21");
            if ((*chromIT).second.count()>0)
                writeBinDensity((*chromIT).second, outFile, binSize);
        }
    }
    else
    {
        for (auto chromIT= tagExp.first(); chromIT!=tagExp.last(); chromIT++)
        {
            //auto tagChrom= tagExp.getpChrom("chr21");
            if ((*chromIT).second.count()>0)
                writeDensityFromFile((*chromIT).second, outFile, bedgraph);
        }
    }
    }
    catch (elem_throw & e)
    {
        const uTags* errorTagPoint;
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
    catch(...)
    {
        cerr << "Caught something else" <<endl;
    }

}


void writeBinDensity( uTagsChrom& tagChrom, std::ostream& out, int binSize)
{

    vector<long int> densityValues;
    densityValues.resize(tagChrom.getChromSize());
    string chrName= tagChrom.getChr();
    //  auto func = ([&] (uTags Elem){        densityValues.at(Elem.getStart());   } );
    /*  const_cast<uTagsChrom*> (&tagChrom)->applyOnAllSites([&] (uTags Elem)
      {
          for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
          {
              densityValues.at(i)++;
          }
      }
    ); */

    for (int j=0; j <((int)densityValues.size()/binSize); j++ )
    {
        int start=j*binSize;
        int end=j*binSize+binSize;
        int tagcount=0;
        tagcount= tagChrom.getSubsetCount(start,end);
                  //  for (int k=start; k< end; k++)
                  // {
                  //  tagcount+=densityValues.at(k);
                  //  }
        out<<chrName << "\t" <<start << "\t" <<end <<"\t" << "a" << "\t" <<tagcount <<endl;
    }

}



void writeDensityFromFile(const uTagsChrom& tagChrom, std::ostream& out , bool bedgraph)
{

    vector<long int> densityValues;
    densityValues.resize(tagChrom.getChromSize());
    //  auto func = ([&] (uTags Elem){        densityValues.at(Elem.getStart());   } );

    const_cast<uTagsChrom*> (&tagChrom)->applyOnAllSites([&] (uTags Elem)
    {
        for (int i=Elem.getStart(); i<Elem.getEnd(); i++)
        {
            // cout << i << endl;
            densityValues.at(i)++;
        }
    }

                                                        );
    auto ourWig=vectorToWig(densityValues,tagChrom.getChr() );
    if (bedgraph)
        writeWigAsBedgraph(ourWig, out);
    else
        writeWig(ourWig, out);

}

void writeWigAsBedgraph(const vector<wigData> & ourData, std::ostream& out)
{

    for(auto& wigValue:ourData )
    {
        out <<wigValue.chr<< "\t" << wigValue.position<< "\t" << (wigValue.position+wigValue.span-1)<< "\t" << wigValue.value<<endl;
    }

}

void writeWig(const vector<wigData> & ourData, std::ostream& out)
{

    const string step="variableStep chrom=";
    const string span="  span=";
    auto curSpan= 0;
    for(auto& wigValue:ourData )
    {
        if (wigValue.span!=curSpan)
        {
            out <<step<<wigValue.chr<<span<<wigValue.span<<endl;
            curSpan=wigValue.span;
        }
        out << wigValue.position << "\t" << wigValue.value << endl;
    }

}

vector<wigData> vectorToWig(vector<long int> densityVector, string chrom)
{

    wigData tempData;
    tempData.position=0;
    tempData.chr=chrom;
    tempData.value=densityVector.at(0);
    vector<wigData> returnVec;
    long long int span=0;
    long long int pos=0;;
    for(auto& densityValue:densityVector )
    {
        if (densityValue!=tempData.value)
        {
            tempData.span=span;
            returnVec.push_back(tempData);
            /**< New data */

            tempData.value = densityValue;
            tempData.position=pos;
            span=1;
        }
        else
            span++;

        /**< Next position */
        pos++;
    }

    if (tempData.span>1)
        returnVec.push_back(tempData);

    return returnVec;
}