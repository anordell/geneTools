#ifndef DENSITYFROMFILE_H_INCLUDED
#define DENSITYFROMFILE_H_INCLUDED

#include <vector>
#include <string>
#include <iostream>

namespace NGS{
    class uTagsChrom;
}
struct wigData
{
    std::string chr;
    long long int position;
    long long int span;
    float value;
};
void densityFromFile(int argc, char* argv[]);
void writeDensityFromFile(const NGS::uTagsChrom& tagChrom, std::ostream& out, bool bedgraph);
void writeWig(const std::vector<wigData> & ourData, std::ostream& out);
std::vector<wigData> vectorToWig(std::vector<long int> densityVector, std::string chrom);
void writeBinDensity(NGS::uTagsChrom& tagChrom, std::ostream& out , int binSize);
void writeWigAsBedgraph(const std::vector<wigData> & ourData, std::ostream& out);
#endif // DENSITYFROMFILE_H_INCLUDED
