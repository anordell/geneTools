#include <bedFile.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{

   bedFile resultBed;

   string input, output;

   input = argv[1];
   output = argv[2];


   // input = "refSeq.bed";
   //output = "happy.txt";

    resultBed.readTags(input);
    resultBed.mapTagstoChrom();
    resultBed.sortChromTags();
    resultBed.outputBedTags(output);

    return 0;
}
