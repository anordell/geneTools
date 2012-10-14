#ifndef GENERATESIGNAL_H_INCLUDED
#define GENERATESIGNAL_H_INCLUDED

#include <vector>
#include <string>

namespace NGS{
    class uRegionExperiment;
    class uTagsExperiment;
}
std::vector<float> binSignal(const std::vector<float> & input_Signal, const int & binSize);
void generateSignal(int argc, char* argv[]);
void setRegionsSize(NGS::uRegionExperiment & ourRegionExp, int binSize);
std::vector<float> getAvgSignal(NGS::uRegionExperiment & ourRegionExp, int binSize);

std::vector<float> normRPM(std::vector<float> input_Signal,const NGS::uTagsExperiment & ourTags);
std::vector<float> getSDSignal(NGS::uRegionExperiment & ourRegionExp, int binSize);

#endif // GENERATESIGNAL_H_INCLUDED
