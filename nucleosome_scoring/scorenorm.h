#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include "uNucleoBin.h"
#include "uRegion.h"
#include "functions.h"

void scoreAndNormalize(int argc, char* argv[]);

void decompose(std::string pathname,std::string chromName, int threshold, bool isComplete, statsStruct pSD);
std::vector<uRegion>  getAndNormRegion(uTagsExperiment* pTags,std::string pChromname,  int pThreshold, bool pIsComplete, statsStruct pSD );
void writeRegNuclAsBedGraph(std::vector<uRegion> * vecRegion, std::string ofRegPath);
