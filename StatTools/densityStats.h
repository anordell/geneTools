#ifndef DENSITYSTATS_H_INCLUDED
#define DENSITYSTATS_H_INCLUDED
#include <iostream>
#include <fstream>
 struct StatsStructure{
     float mean;
     float sd;
     float q1;
     float median;
     float q3;
 };


class uTagsExperiment;
class uRegionExperiment;


StatsStructure getCoverageDensityStats(uRegionExperiment& regionExp, uTagsExperiment& tagExp);
StatsStructure getNormDensityStats(uRegionExperiment& regionExp);
StatsStructure getDensityStats(uRegionExperiment& regionExp);
void densityStats(int argc, char* argv[]);

void writeData(StatsStructure statResult, std::ostream & out);

#endif // DENSITYSTATS_H_INCLUDED
