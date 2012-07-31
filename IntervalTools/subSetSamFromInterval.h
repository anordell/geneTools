#ifndef _SUBSETSAMFROMINTERVAL_H_INCLUDED
#define _SUBSETSAMFROMINTERVAL_H_INCLUDED
#include <string>

class uTagsExperiment;
class uRegionExperiment;

void subSetSamFromInterval(int argc, char **argv);
uTagsExperiment loadSamData(std::string path);

uRegionExperiment loadRegionData(std::string path);
#endif // _SUBSETSAMFROMINTERVAL_H_INCLUDED
