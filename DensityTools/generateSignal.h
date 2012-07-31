#ifndef GENERATESIGNAL_H_INCLUDED
#define GENERATESIGNAL_H_INCLUDED

#include <vector>
#include <string>


class uRegionExperiment;
class uTagsExperiment;

std::vector<float> binSignal(const std::vector<float> & input_Signal, const int & binSize);
void generateSignal(int argc, char* argv[]);
void setRegionsSize(uRegionExperiment & ourRegionExp, int binSize);
std::vector<float> getAvgSignal(uRegionExperiment & ourRegionExp, int binSize);

std::vector<float> normRPM(std::vector<float> input_Signal,const uTagsExperiment & ourTags);
std::vector<float> getSDSignal(uRegionExperiment & ourRegionExp, int binSize);

#endif // GENERATESIGNAL_H_INCLUDED
