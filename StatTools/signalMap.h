#ifndef SIGNALMAP_H_INCLUDED
#define SIGNALMAP_H_INCLUDED
#include <string>
#include <vector>
#include <map>
#include <iostream>

void signalMap(int argc, char* argv[]);
std::map<std::string,std::vector<float>> loadSignalDataFromFile(std::ifstream & readfile);
void writeHeatMaptoPS(const std::vector<std::vector<float>> & heatMap, const std::string & psName,const std::vector<std::string> & legendVector, const std::string & titlename = "null" );
std::vector<std::vector<float>> generateDistMatrice(std::vector<std::vector<float>> listOfSignals);
#endif // SIGNALMAP_H_INCLUDED
