#ifndef ALPHADISTRIBUTION_H_INCLUDED
#define ALPHADISTRIBUTION_H_INCLUDED

#include <map>
#include <vector>
#include <string>
#include <ostream>

#define AlphaMap std::map<std::string, std::map<char,int>>

void AlphaDistribution(int argc, char* argv[]);
AlphaMap loadData(std::ifstream  &inputStream);
void processMap(const AlphaMap & alphaData, std::vector<char> ignoreVecChar);
#endif // ALPHADISTRIBUTION_H_INCLUDED
