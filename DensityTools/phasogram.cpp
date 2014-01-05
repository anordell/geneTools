#include "phasogram.h"

#include <string>
#include <string.h>
#include <fstream>
#include "uTags.h"
#include "utility/utility.h"
#include "uRegion.h"
using namespace NGS;
using namespace std;
/** \brief Taking a tab file and Sam file as input, we return the number of tags from the sam file overlapping our tab regions.
 *
 * \param argc int
 * \param argv[] char*
 * \return void
 *
 */
void parsePhasogram(int argc, char* argv[])
{
    int graphSize=0, pileSize=0;
    string pathname="";
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-f")==0)
        {
            pathname = argv[i + 1];
        }
        else if (strcmp(argv[i],"-s")==0)
        {
            graphSize = atoi(argv[i + 1]);
        }
        else if (strcmp(argv[i],"-p")==0)
        {
            pileSize = atoi(argv[i + 1]);
        }
    }
    //If we did not input file or chrom name
    if ((pathname.size()==0)||(graphSize==0)||(pileSize==0))
    {
        cerr<<"Program signature is -f <filepath> -s <Graph Size> -p <Pile Size>";
        return;
    }

    ifstream inputStream(pathname);
    if( !(inputStream))
    {
        cerr <<"Invalid file name for Sam File. File name is "<< pathname <<endl;
        return;
    }

    uTagsExperiment curExperiment;
    uParser samParser(&inputStream,"SAM");
    curExperiment.loadWithParser(samParser);
    auto result=phasogram(curExperiment,pileSize, graphSize);

    for (int x:result)
    {
        cout << x << " ";
    }
    return;
}


vector<int> phasogram(const uTagsExperiment & tagExp, int pileSize, int graphSize)
{

//    tagExp.applyOnAllChroms()
    const uTagsExperiment* b=&tagExp;

    uTagsExperiment* notConst= const_cast<uTagsExperiment* >(b);

    // auto pChrom =  notConst->getpChrom("chr18");
    vector<int> phasogram;
    phasogram.resize(graphSize);
    for (auto chromIT= notConst->begin(); chromIT!=notConst->end(); chromIT++)
    {
        //auto tagChrom= tagExp.getpChrom("chr21");
        if ((*chromIT).second.count()>0)
        {
            auto result = mapTagNGStoDensity((*chromIT).second);
            cerr << "Entering Phasogram " <<endl;
            for(int k=0; k<(long long int)result.size(); k++ )
            {
                if (k%300000==0)
                    cerr << "done pos " << k << endl;
                for(int i=1; i<graphSize; i++)
                {
                    if ((k+i)<(int)result.size())
                        if (result.at(k)>pileSize && result.at(k+i) > pileSize)
                        {
                            phasogram.at(i)++;
                            break;
                        }

                }

            }
        }
    }
    return phasogram;
}

vector<int> mapTagNGStoDensity( uTagsChrom & tagChrom)
{
    vector<int> returnVec;
    if (( (tagChrom.getChromSize()+1)*2) > (long long int)returnVec.max_size() )
    {
        cerr<< "Not enough memory for vector allocation";
        abort();
    }
    //We resize multiple times to avoid blowing our stack?
    returnVec.resize(tagChrom.getChromSize()+1);
    uTags tempTag;
    int count=0;
    bool end=false;

    if (tagChrom.count()>count)
        tempTag=tagChrom.getSite(count);
    else
        end =true;
    while(!end)
    {
        count++;

        int pos=0;
        //Is strand is positive ( should be always for completed fragments ), get central position based on paired lenght or fragment lenght
        //If not, do nothing
        if (tempTag.getStrand()==StrandDir::FORWARD)
        {
            pos= (tempTag.getStart());
            //If central position is beyond our registered area, error. This happens as BWA something maps over our reference?
            if (pos >=(long long int)returnVec.size() )
            {
                cerr <<"Invalid tag position at pos " << pos <<endl;
            }
            else //Density ++ at central point
                returnVec.at(pos)=returnVec.at(pos)+1;
        }

        //Since a Vector is zero based, this works! Syntax a bit unclear thought

        if (tagChrom.count()>count)
            tempTag=tagChrom.getSite(count);
        else
            end =true;
    }
    //Random info output
    cerr << "Loaded " << count << " tags info" << endl;
    return returnVec;
}


