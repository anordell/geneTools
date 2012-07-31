#ifndef EXTEND_INCLUDED
#define EXTEND_INCLUDED
#include <string>



enum class ExtendType
{
    LEFT,RIGHT,BI
};

class uIntervalExperiment;

void extend(int argc, char **argv);
uIntervalExperiment loadData(std::string path);
uIntervalExperiment extendData(uIntervalExperiment intervalData, ExtendType parType, unsigned long int pShiftSize);
#endif // EXTEND_INCLUDED
