
#ifndef SCOREH_INCLUDED
#define SCORE_H_INCLUDED
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include "uNucleoBin.h"

void score(int argc, char* argv[]);
void mainScoring(std::string pathname, uNucleoBin &ourBin, int threshold, bool isComplete);

#endif
