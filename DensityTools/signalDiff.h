#ifndef SIGNALDIFF_H_INCLUDED
#define SIGNALDIFF_H_INCLUDED

namespace NGS{

class uRegionExperiment;
}

void generateDistanceScores( NGS::uRegionExperiment & regionExpA, NGS::uRegionExperiment & regionExpB);
void signalDiff(int argc, char* argv[]);

#endif // SIGNALDIFF_H_INCLUDED
