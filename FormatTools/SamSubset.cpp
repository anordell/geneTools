

#include "SamSubset.h"
//This writes, in Sam format, the tags from a subset of our data.
using namespace std;
void getSamSubset(uRegion ourRegion, uTagsExperiment* ourExp, std::ostream& output){

    uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>,uGenericNGS> testExp;
    auto yar = (ourExp->getOverlapping(ourRegion.getChr(),ourRegion.getStart(),ourRegion.getEnd()));
   // tempExp=static_cast<uTagsExperiment>(ourExp->getOverlapping(ourRegion.getChr(),ourRegion.getStart(),ourRegion.getEnd()));
    auto tempExp =static_cast<uTagsExperiment*>(&yar);

    tempExp->writeSamToOutput(output);
   // tempExp.writeSamToOutput(output);
}
