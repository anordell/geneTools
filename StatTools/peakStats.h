#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <time.h>
#include <thread>
#include <random>
#include <string.h>


enum class GenomicFileType;


void peakStatistics(int argc, char* argv[]);
void loadData(std::ifstream&, const GenomicFileType,std::ostream& out);
