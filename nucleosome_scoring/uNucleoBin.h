#ifndef UNUCLEOBIN_H_INCLUDED
#define UNUCLEOBIN_H_INCLUDED
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "uNucleoBin.h"
#include "uRegion.h"
#include "functions.h"
//#include <utilityh>
struct statsStruct;
class uNucleoBin{

    public:

      void generateStartDensityBin(uTagsExperiment* ourTagExp, int chromSize, int threshold, bool ignoreEmpty);
      void generateDensityBin(uTagsExperiment* ourTagExp, int chromSize, int threshold, bool isComplete, bool ignoreEmpty=false);
      void generateDensityBinFromComplete(uTagsExperiment* ourTagExp, int chromSize, int threshold, bool ignoreEmpty=false);
      void writeDensitytoFile(std::ostream &out);
      void writebedGraph(std::ostream &out, float threshold=0);
      void loadScoreFromBedGraph(std::ifstream &in);

      std::string getChromName(){return chromName;}
      void setChromName(std::string chrom){chromName=chrom;}
      void generateSMap();
      void generateDMap();


      std::vector<int> * getpMapCount(){return &DensityMapCount;};
      std::vector<int>::iterator CountBegin(){return DensityMapCount.begin();};
      std::vector<int>::iterator CountEnd(){return DensityMapCount.end();};


      std::vector<float> * getpMapScore(){return &DensityMapScore;};
      std::vector<float>::iterator Scorebegin(){return DensityMapScore.begin();};
      std::vector<float>::iterator Scoreend(){return DensityMapScore.end();};
      float getScore(int pos){return DensityMapScore.at(pos);};
      void writeCounttoBedGraph(std::ostream &out, int range);
      void GaussianNormSMap(float sd, float mean);
      uNucleoBin(){};
      statsStruct getSD();
       uNucleoBin(uNucleoBin&& other)
        : DensityMapCount( std::move(other.DensityMapCount) )
        , DensityMapScore( std::move(other.DensityMapScore) )
        , chromName( other.chromName)
        , m_ignoreThreshold(other.m_ignoreThreshold)
        , m_ignoreEmpty(other.m_ignoreEmpty)
        , L(other.L)
    {}

       uNucleoBin(const uNucleoBin &other)
        : DensityMapCount( other.DensityMapCount )
        , DensityMapScore( other.DensityMapScore )
        , chromName( other.chromName)
        , m_ignoreThreshold(other.m_ignoreThreshold)
        , m_ignoreEmpty(other.m_ignoreEmpty)
        , L(other.L)
    {}

    uNucleoBin& operator=(const uNucleoBin& other)
    {
        DensityMapCount=other.DensityMapCount;
        DensityMapScore=other.DensityMapScore;
        chromName=other.chromName;
        m_ignoreThreshold=other.m_ignoreThreshold;
        m_ignoreEmpty=other.m_ignoreEmpty;
        L=other.L;

       return *this;
    }

    uNucleoBin& operator=(uNucleoBin&& other)
    {
        DensityMapCount= std::move(other.DensityMapCount);
        DensityMapScore= std::move(other.DensityMapScore);
        chromName= std::move(other.chromName);
        m_ignoreThreshold=other.m_ignoreThreshold;
        m_ignoreEmpty=other.m_ignoreEmpty;
        L=other.L;
       return *this;
    }

       std::vector<int> DensityMapCount;
        std::vector<float> DensityMapScore;
    //  vector<bedScores> getbedVec();
    private:

        struct startCount{
            int count;
            float score;
        };

        //library size
        static const int l=150;
        //Start Count


        std::string chromName;

        int m_ignoreThreshold;
        bool m_ignoreEmpty;
        //ChromSize;
        int L;
        //number of dyads at position j
        int d(int j);

        //Smoothing Kernel
        float  K(float u, float w);

        //Kernel-smoothed Dyad count
        float D(int i, int w);

        //Stringency at position i
        float S(int i, int w);

        static const int w=30;

};


namespace clustering
{

    std::vector<float> hauftsmanTwoBins(const uNucleoBin &  binA,const uNucleoBin & binB );

}


#endif // UNUCLEOBIN_H_INCLUDED
