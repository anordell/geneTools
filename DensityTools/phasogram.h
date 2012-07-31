#ifndef PHASOGRAM_INCLUDED
#define PHASOGRAM_INCLUDED

#include "uTags.h"
#include <vector>
#include <string>

class uTagsExperiment;


std::vector<int> phasogram(const uTagsExperiment & tagExp, int pileSize, int graphSize);
std::vector<int> mapTagNGStoDensity( uTagsChrom & tagChrom);
void parsePhasogram(int argc, char* argv[]);


#endif // PHASOGRAM_INCLUDED
