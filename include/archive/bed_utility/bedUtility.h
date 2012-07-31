#ifndef BEDUTILITY_H_INCLUDED
#define BEDUTILITY_H_INCLUDED

#define RETURN_PARTIAL 0
#define RETURN_COMPLETE 1
#define RETURN_COMPLETE_FLANK 2

struct barrierResults{
std::vector<TagRead> barrierTags[25];
std::vector<TagRead> centerTags[25];
std::vector<TagRead> otherTags[25];
};


bedFile overlapBedfiles(bedFile* fileA, bedFile* fileB, int overlapcount,int distance, int overlap_para);

barrierResults barrierOverlap(bedFile inputData, bedFile comparisonData, int distance);

bedFile substractBedfiles(bedFile* fileA, bedFile* fileB, int overlapcount,int distance);

TagRead isoverlapTag(TagRead tagA, TagRead tagB, int overlapcount, int distance,int overlap_para);

     //Function to swap from internal to string representation of a chrom number ( 0 = Chrom 1, etc )
     std::string getChromfromint(int chromNum);
     int getChromfromString(std::string);


#endif // BEDUTILITY_H_INCLUDED
