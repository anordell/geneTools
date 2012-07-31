#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "uNucleoBin.h"
#include "uRegion.h"
#include "uFunctions.h"

struct bedScores{
std::string chr;
int position;
int end;
float score;
};

struct statsStruct{
float sd;
int sum;
float mean;
};
/**< Forward dec */
class uNucleoBin;
 /**< End forward */
float normScore(float sd, float mean);
uTagsExperiment decompPass(uTagsChrom* ourTagChrom, std::vector<uRegion> vecRegions);
void decomp(int arg, char* argv[]);
std::vector<uRegion> filterNucleoMaxima(uTagsExperiment* ourTagExp, std::string chr, std::vector<bedScores> &vecScores);
std::vector<bedScores> loadbedGraph(std::ifstream& inputStream);
std::vector<bedScores> loadbedGraph(std::string bedpath);
statsStruct returnSdfromRegions(std::vector<uRegion> regionVec);
std::vector<uRegion> bedgraphToRegions(uTagsExperiment ourExp,const std::vector<bedScores> vecBedscore, int size, std::string chromName);
std::vector<bedScores> returnMaximaBedVec(const std::vector<bedScores> bedVec);
std::vector<uRegion> makeRegionsfromBin(uNucleoBin& ourBin,uTagsExperiment& ourExp,  float threshold=0);
std::vector<bedScores> getMaximaFromBedscores(const std::vector<bedScores> ourBedScores);
void Normalize(float sd, float mean,std::vector<uRegion> & ourRegions);
std::vector<uRegion> filterRegions(std::vector<uRegion> ourRegions, float threshold);
std::vector<bedScores> getBedFromEnlargedRegions(std::vector<uRegion> & ourRegions, int enlargeSize);
uNucleoBin loadDensityFromTags(uTagsExperiment& ourExp, int threshold,std::string chromname, bool isComplete);


#endif
