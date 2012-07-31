#ifndef DIFFPLOT_H_INCLUDED
#define DIFFPLOT_H_INCLUDED
#include <vector>
#include "uRegion.h"


//class uRegionExperiment;
class uTagsExperiment;
//class uRegion;

struct Elem_score{
float score;
uRegion Elem;
};

std::vector<Elem_score> compareSignals(uRegionExperiment & markA,uRegionExperiment & markB);
void diffPlot(int argc, char* argv[]);
void setRegionsSize(uRegionExperiment & ourRegionExp, int extendSize);

#endif // DIFFPLOT_H_INCLUDED
