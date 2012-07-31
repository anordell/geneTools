#include "uNucleoBin.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include "functions.h"
#include "clustering.h"
using namespace std;
//This should be deprecated normally!
void uNucleoBin::generateDensityBinFromComplete(uTagsExperiment* ourTagExp, int chromSize, int threshold, bool ignoreEmpty)
{
    int pos;

    //  if (ourTagExp.)

    m_ignoreThreshold  =threshold;
    m_ignoreEmpty= ignoreEmpty;

    if (( (chromSize+1)*2) > (long long int)DensityMapCount.max_size() )
    {
        cerr<< "Not enough memory for vector allocation";
        abort();
    }
    cerr << "Requesting vector(twice) resize to " << (chromSize+1) << " int "<<endl;
    cerr << "Max size should be "  << DensityMapCount.max_size() << endl;

    //We resize multiple times to avoid blowing our stack?
    DensityMapCount.resize((chromSize+1));
    DensityMapScore.resize((chromSize+1));
    cerr << "Allocated memory" << endl;

    L= chromSize;
    uTags tempTag;
    int count=0;
    while(ourTagExp->isEndfile()==false )
    {
        count++;
        tempTag= ourTagExp->nextSamLine();
        pos=0;
        if (tempTag.getChr()==this->getChromName())
        {
            if (tempTag.getStrand()=='+')
            {
                pos= ( tempTag.getStart()+ (tempTag.getLenght()/2) );
            }
        }
        if (pos >=(int)DensityMapCount.size() )
        {
            cerr <<"Invalid tag at position" << pos <<endl;
            cerr <<"Max Size is" <<(int)DensityMapCount.size() <<endl;
            cerr << tempTag.getStart() <<endl;
            cerr << tempTag.getEnd() <<endl;
            cerr << tempTag.getStrand() <<endl;
            cerr << "Is PE " << tempTag.isPE()<<endl;
            cerr << tempTag.getPeLenght()/2 <<endl;
        }
        else
        {
            if(pos)
                DensityMapCount.at(pos)=DensityMapCount.at(pos)+1;
        }
        if (count%200000==0)
        {
            cerr <<"Processed " << count  << " tags "<<endl;
        }
    }
    //Generate our Stringency map
    cerr << "Loaded " << count << " tags info" << endl;
}

void uNucleoBin::generateStartDensityBin(uTagsExperiment* ourTagExp, int chromSize, int threshold, bool ignoreEmpty)
{
    int pos;

    m_ignoreThreshold  =threshold;
    m_ignoreEmpty= ignoreEmpty;

    if (( (chromSize+1)*2) > (long long int)DensityMapCount.max_size() )
    {
        cerr<< "Not enough memory for vector allocation";
        abort();
    }
    cerr << "Requesting vector(twice) resize to " << (chromSize+1) << " int "<<endl;
    cerr << "Max size should be "  << DensityMapCount.max_size() << endl;

    //We resize multiple times to avoid blowing our stack?
    DensityMapCount.resize((chromSize+1));
    DensityMapScore.resize((chromSize+1));
    cerr << "Allocated memory" << endl;
    L= chromSize;
    uTags tempTag;
    int count=0;
    bool end=false;
    uTagsChrom* pChrom =nullptr;

    if(ourTagExp->isModeGradual())
    {

        if (ourTagExp->isEndfile())
            end=true;
        else
            tempTag= ourTagExp->nextSamLine();
    }
    else
    {
        pChrom=ourTagExp->getpChrom(this->getChromName());
        if (pChrom->count()>count)
            tempTag=pChrom->getSite(count);
        else
            end =true;
    }
    while(!end)
    {
        count++;

        if (tempTag.getChr()==this->getChromName())
        {
            pos=0;
            //Is strand is positive ( should be always for completed fragments ), get central position based on paired lenght or fragment lenght
            //If not, do nothing
            if (tempTag.getStrand()=='+')
            {
                pos= ( tempTag.getStart());

                //If central position is beyond our registered area, error. This happens as BWA something maps over our reference?
                if (pos >=(int)DensityMapCount.size() )
                {
                    string a;
                    cerr <<"Invalid tag position at pos " << pos <<endl;
                    cerr << tempTag.getName() <<endl;
                    cerr << tempTag.getChr() <<endl;
                    cerr << tempTag.getStart() <<endl;
                    cerr << tempTag.getEnd() <<endl;
                    cerr << tempTag.getStrand() <<endl;
                    cerr << "Is PE " << tempTag.isPE()<<endl;
                    cerr << tempTag.getPeLenght()/2 <<endl;

                }
                else //Density ++ at central point
                    DensityMapCount.at(pos)=DensityMapCount.at(pos)+1;
            }
        }
        if (count%200000==0)
        {
            cerr <<"Processed " << count  << " tags "<<endl;
        }

        //Load our next tag
        if(ourTagExp->isModeGradual())
        {
            if (ourTagExp->isEndfile())
                end=true;
            else
                tempTag= ourTagExp->nextSamLine();
        }
        else
        {
            //Since a Vector is zero based, this works! Syntax a bit unclear thought
            if (pChrom)
            {
                if (pChrom->count()>count)
                    tempTag=pChrom->getSite(count);
                else
                    end =true;
            }
        }
    }
    //Random info output
    cerr << "Loaded " << count << " tags info" << endl;
    cerr << "Density map is of size " <<DensityMapCount.size() <<endl;
}

/** \brief Count every tag that overlaps every position and store
 *
 * \param ourTagExp uTagsExperiment* The tags to use
 * \param chromSize int Size of the chromosome ( should be implicit? )
 * \param threshold int Below this count, we ignore
 * \param isComplete bool False is working with Paired end data, true if we have already completed our PE data
 * \param ignoreEmpty bool Parameter for the future???
 * \return void
 *
 */
void uNucleoBin::generateDensityBin(uTagsExperiment* ourTagExp, int chromSize, int threshold, bool isComplete, bool ignoreEmpty)
{
    long long int pos;

    m_ignoreThreshold  =threshold;
    m_ignoreEmpty= ignoreEmpty;

    if (( (chromSize+1)*2) > (long long int)DensityMapCount.max_size() )
    {
        cerr<< "Not enough memory for vector allocation";
        abort();
    }
    cerr << "Requesting vector(twice) resize to " << (chromSize+1) << " int "<<endl;
    cerr << "Max size should be "  << DensityMapCount.max_size() << endl;

    //We resize multiple times to avoid blowing our stack?
    DensityMapCount.resize((chromSize+1));
    DensityMapScore.resize((chromSize+1));
    cerr << "Allocated memory" << endl;
    L= chromSize;
    uTags tempTag;
    int count=0;
    bool end=false;
    uTagsChrom* pChrom=nullptr;
    ofstream temp("TagErrors.txt");
    if(ourTagExp->isModeGradual())
    {

        if (ourTagExp->isEndfile())
            end=true;
        else
            tempTag= ourTagExp->nextSamLine();
    }
    else
    {
        pChrom=ourTagExp->getpChrom(this->getChromName());
        if (pChrom->count()>count)
            tempTag=pChrom->getSite(count);
        else
            end =true;
    }
    while(!end)
    {
        count++;

        if (tempTag.getChr()==this->getChromName())
        {
            pos=0;
            //Is strand is positive ( should be always for completed fragments ), get central position based on paired lenght or fragment lenght
            //If not, do nothing

            if (tempTag.getStrand()=='+')
            {
                if (isComplete)
                    pos= ( tempTag.getStart()+ (tempTag.getLenght()/2) );
                else
                    if (tempTag.isPE())
                        pos= ( tempTag.getStart()+ (tempTag.getPeLenght()/2) );
                    else
                        pos=-1;

                //If central position is beyond our registered area, error. This happens as BWA something maps over our reference?
                if ((pos >=(int)DensityMapCount.size() )||(pos<0))
                {
                    string a;
                    temp <<"Invalid tag position at pos " << pos <<endl;
                    temp << tempTag.getName() <<endl;
                    temp << tempTag.getChr() <<endl;
                    temp << tempTag.getStart() <<endl;
                    temp << tempTag.getEnd() <<endl;
                    temp << tempTag.getStrand() <<endl;
                    temp << "Is PE " << tempTag.isPE()<<endl;
                    temp << tempTag.getPeLenght()/2 <<endl;
                }
                else //Density ++ at central point
                {

                    DensityMapCount.at(pos)=DensityMapCount.at(pos)+1;
                }
            }
        }
        if (count%200000==0)
        {
            cerr <<"Processed " << count  << " tags "<<endl;
        }

        //Load our next tag
        if(ourTagExp->isModeGradual())
        {
            if (ourTagExp->isEndfile())
                end=true;
            else{
                tempTag= ourTagExp->nextSamLine();
                if ( (tempTag.getStart()==tempTag.getEnd())&&(tempTag.getStart()==0))
                    end=true;
            }
          //  cerr << "end load" <<endl;
        }
        else
        {
            //Since a Vector is zero based, this works! Syntax a bit unclear thought
            if (pChrom->count()>count)
                tempTag=pChrom->getSite(count);
            else
                end =true;
        }


    }


  //  for (int x=0; x< DensityMapCount.size();x++)
   //    temp << x<< " "<< DensityMapCount.at(x) <<endl;
    //Random info output
    cerr << "Loaded " << count << " tags info" << endl;
    cerr << "Density map is of size " <<DensityMapCount.size() <<endl;
}



/**< Write the scoring map passe nature formale? */
void uNucleoBin::generateSMap()
{
    int i;

        #pragma omp parallel for private(i)
        for(i=1; i<= L; i++ )
        {
            DensityMapScore.at(i)= S(i, w);
        }
        //ofstream score "score.txt";

    cerr << "Finished processing positions"<<endl;

}
//Get d(j)
int uNucleoBin::d(int j)
{
    int value;

    value = DensityMapCount.at(j);

    return value;
}

float uNucleoBin::K(float u, float w)
{

    float fP, sP, tP;
    // I{abs(u)<w} = 1
    if (abs(u)>=w)
        return 0;

    fP= pow((u/w),2);
    sP=(1- fP);
    tP= pow(sP,3);
    return tP;
}
//Kernel Smoothed Dyad count
//We can optimise this.
float uNucleoBin::D(int i, int w)
{

    float sum;
    float count=0;
    sum=0;
    int bot, top, botCheck, topCheck;
    bot = (i-(w*2));
    if ( bot<0)
        bot = 0;

    top =(i+(w*2));
    if (top > L)
        top =(L-1);

    botCheck = (i-(w*5));
    if ( botCheck<0)
        botCheck = 0;

    topCheck =(i+(w*5));
    if (topCheck > L)
        topCheck =(L-1);
    if (m_ignoreThreshold)
    {
        for(int j=botCheck; j<=topCheck; j++)
        {
            count+=d(j);
        }

        if (count <m_ignoreThreshold )
            return 0;
    }

    for(int j=bot; j<=top; j++)
    {
        if (d(j)!=0)
            sum+=( K((i-j), w)*d(j));
    }


    return sum;
}

float uNucleoBin::S(int i, int w)
{

    float top, bottom, B,T;
    //Numerator

    if (m_ignoreEmpty)
    {
        if ((i>=L)||((d(i)==0)))
            return 0;
    }
    top= D(i,w);
    bottom =0;
    float result;
    //If numerator is empty, no need to do this.
    if (top!=0)
    {

        //Denominator
        //Limit the formulae demands.
        B=(i-150);
        if (B<0)
            B=0;

        T=(i+150);
        if (T>L)
            T=L;
        for(int j=B; j<=T; j++)
        {

            bottom+=( (1.09/w)*D(j,w) );
        }
    }

    if (bottom==0)
        result=0;
    else
        result = (top/bottom);

    return result;
}


void uNucleoBin::writeDensitytoFile(std::ostream &out)
{

    cerr <<" Writing data to file";
    for(unsigned int k=0; k< DensityMapCount.size(); k++)
    {
        if (DensityMapScore.at(k)!=0)
        {
            out<< k << "\t" << DensityMapScore.at(k) <<"\n";

        }
    }

}

void uNucleoBin::writeCounttoBedGraph(std::ostream &out, int range)
{
//Header
    out << " type=bedGraph" << endl;
    unsigned int bot, top;
    int count;
    for(unsigned int k=0; k< DensityMapCount.size(); k++)
    {

        bot = k-range;
        top = k+range;
        if (top>=DensityMapCount.size())
            top=DensityMapCount.size();
        if (bot<0)
            bot=0;
        for (unsigned int i=bot; i<= (top); i++ )
            count+=DensityMapCount.at(i);

        if (count!=0)
        {
            out  << this->getChromName() << "\t" << (k-1) << "\t" << k << "\t" << count <<"\n";
        }
        count=0;
    }

}

void uNucleoBin::writebedGraph(std::ostream &out, float threshold)
{

//Header
    out << " type=bedGraph" << endl;

    for(unsigned int k=0; k< DensityMapCount.size(); k++)
    {
        if (DensityMapScore.at(k)>threshold)
        {
            out  << this->getChromName() << "\t" << (k-1) << "\t" << k << "\t" << DensityMapScore.at(k) <<"\n";
        }
    }
}

/** \brief Load our score form an existing file

 *
 * \param in std::istream&
 * \return void
 *
 */
void uNucleoBin::loadScoreFromBedGraph(std::ifstream &in)
{
   auto result = loadbedGraph(in);
   for( bedScores& score : result)
        {
            if ((long long int)DensityMapScore.size()<=score.position)
                DensityMapScore.resize(score.position+1);
            DensityMapScore.at(score.position)=score.score;
        }
}

/**< Functions  */

namespace clustering{



 /** \brief Distance between two equaled sized Score vectors
  *
  * \param binA const uNucleoBin First bim
  * \param binB const uNucleoBin Second bin
  * \return vector<float> List of scores for each bin
  *
  */
 vector<float> hauftsmanTwoBins(const uNucleoBin &  binA,const uNucleoBin & binB )
 {
    vector<float> returnDistance;
    try{
        if (binA.DensityMapScore.size()!=binB.DensityMapScore.size())
            throw 10;

        //For each bin of size Y, measure distance

        for(int i=0; i<(int)binA.DensityMapScore.size()/100; i++)
        {
            vector<float> subsetA,subsetB;

            for(int k=(i*100); k<(i*100+100); k++){
                subsetA.push_back(binA.DensityMapScore.at(k));
                subsetB.push_back(binB.DensityMapScore.at(k));
            }
            auto result= clustering::hauftsman(subsetA, subsetB);
           if (result > 0)
            cout << "chr18" << "\t" <<i*100<< "\t" << (i*100+99) << "\t" << result << endl;
            returnDistance.push_back(result);

        }


     }
     catch(...){
        cerr << "in hauftsmanTwoBins, bins as not same size" <<endl;
     throw;
     }

return returnDistance;
 }



}
/**< Measure gaussian distance and ponder the distance */
void uNucleoBin::GaussianNormSMap(float sd, float mean)
{

    //We have our mean and density deviation
    float mult;

  //  DensityMapScore.at(k)
    int i=0;
    //Ponder every score based on it's density. Shoudl we be checking empty??
    for ( float& curscore :DensityMapScore )
    {

        if (DensityMapCount.at(i)<(mean))
        {
            curscore=((utility::gaussianSim(DensityMapCount.at(i),(mean),(sd/4)))*curscore);
        }
        else
        {
            mult= (1-utility::gaussianSim(DensityMapCount.at(i),(mean),(sd/4)));
            curscore=(1+mult)*curscore;
        }
        i++;
    }

}


/**<  Get SD and mean of our density data*/
statsStruct uNucleoBin::getSD()
{
    vector<bedScores>::iterator it;
    vector<float>::iterator floatit;
    vector<float>::iterator intit;
    vector<uRegion>::iterator regit,ic;
    float sum=0, mean=0, sd=0,sumsq=0;

/**< Listen. This is not sd of our data. Is is the sd of the nucleosomeal density for the region.l
    So for each region that has at least 1 tag, we sum over 60 bp the density of tags and make an entry with that.
    This means, they will often overlap (sum of 0-60 then sum 1-61 ), etc
    After, we have a vector with each position being the density of one 60 bp region.
    We measure our SD on THAT!
 */

    /**< For every point.. make a potential nucleosome *///Establish our score passed on density average and sigma
    vector<float> deviations;
    vector<float> nucleoSum;
    int k=0;

    for (float score: DensityMapScore )
    {
       /**< 30 before and after */
        if (score >0)
        {
           int bot = k-30;
           int top= k+30;
           if (bot<0)
                bot =0;
            if (top >L)
                top=(L-1);
            int tempsum=0;
            for (long long int z=bot; z<=top; z++){

                  tempsum += DensityMapCount.at(z);
            }
            sum +=tempsum;
             nucleoSum.push_back(tempsum);
        }
        k++;
    }
    mean = sum/(nucleoSum.size());
    for (int density: nucleoSum )
    {
        deviations.push_back(density-mean);

    }

    for (intit=deviations.begin(); intit < deviations.end(); intit++ )
    {
        *intit=(pow(*intit,2));
        sumsq+=*intit;
    }

    sd = sumsq/(deviations.size()-1);
    sd= sqrt(sd);


    cerr << "Region count is " << nucleoSum.size() <<endl;
    cerr <<"sum  is"<< sum <<endl;
    cerr <<"sumQ  is"<< sumsq <<endl;
    cerr << "deviation size is" << (deviations.size()-1) <<endl;


    statsStruct returnValue;

    returnValue.sd=sd;
    returnValue.sum=sumsq;
    returnValue.mean=mean;

    return returnValue;
}
