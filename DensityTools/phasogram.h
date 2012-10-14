#ifndef PHASOGRAM_INCLUDED
#define PHASOGRAM_INCLUDED

#include "uTags.h"
#include <vector>
#include <string>

//class NGS::uTagsExperiment;

std::vector<int> phasogram(const NGS::uTagsExperiment & tagExp, int pileSize, int graphSize);
std::vector<int> mapTagNGStoDensity( NGS::uTagsChrom & tagChrom);
void parsePhasogram(int argc, char* argv[]);


#endif // PHASOGRAM_INCLUDED
