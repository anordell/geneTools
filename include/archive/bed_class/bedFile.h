#ifndef _BEDFILE
#define _BEDFILE_
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "utility.h"
#include "chromConstants.h"

struct TagIsland{

 int iChrom;
 //Number of tags
 int iCount;
 int iStart;
 int iStop;
 int iWindowCount;
};


struct TagRead{

 //Use for tags and or other stuff
 int chrom;
 int chromStart;
 int chromEnd;

 std::string geneID;
};


struct ChromWindow{

 //Chromosome
 int wStart;
 int wCount; //Number of aligned tags in the window.
 int wSize; //We keep for irregular steps.
};

struct ChromSum{
    pairestd::vector<TagRead> tList; //Tags of the chromosome
    std::vector<ChromWindow> wList; //Windows of the Chromosome
    std::vector<TagIsland>  iList;//Islands on the Chromose;
    int ChromNum; //Chromosome number -1 ( so Chr1= 0 ). See Constant def for CHROMM, X and Y
    int ChromSize;

    //Default size of ours windows, -1 if irregulars sizes.
    int wSize;

};

typedef std::vector<TagRead> Vec_Tag;


class bedFile{


  private:
    //Vector of all tags we want to use/map.
    Vec_Tag tagfile;
    //Have we mapped
     bool isTagMapped;
     bool isWindowMapped;
     bool isIslandMapped;
    //Summary of information per Chromosomes
     ChromSum genomeSumArray[25];


    //Functions
     bool initiWindowSpace(int wSize);
     //Function to swap from internal to string representation of a chrom number ( 0 = Chrom 1, etc )
     std::string getChromfromint(int chromNum);
     int getChromfromString(std::string);
     //Swaps two tags in a tag list.
     void swapTags(Vec_Tag* tag_data,int posA, int posB);


    //Sort functions
     int partition(ChromSum* chromSum,  int &pivot,  int startIndex,  int endIndex);
     void QuickSortChrom(ChromSum* pCurChrom,  int startIndex,  int endIndex);

  public:

  bedFile();
  ~bedFile();
     //Asks for file name
      bool openReadTags();
     //Does not
      bool readTags(std::string);
      //Map tags into windows?
      bool mapTagstoWindows(int wSize) ;
      //Map tags to there chromosomes, no windows or islands
      bool mapTagstoChrom();
      //After mapping into windows, call this to create islands.
      bool countIslands(const int minislanddensity, const int emptywindowLimit, const int minwindowdensity);

      //Functions to manually add data.

      void addRead(TagRead);
      void addReadtoChrom(TagRead, int);

    //    Sort our loaded chromosome info by start.
      void sortChromTags();



        //Output functions
      void outputBedTags(std::string);
      void outputWindows(std::string o_name);
      void outputIslands(std::string);

      int tagCount(){return tagfile.size();}

      void reset();

      ChromSum* getChrom(int chromNumber);
      const Vec_Tag* getTagVec();

};


#endif
