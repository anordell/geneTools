#include <bedFile.h>
#include <bedUtility.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;


int main(int argc, char* argv[])
{
 bedFile bedA,bedB, resultBed;

int overlapcount, distance;

string fileA, fileB, fileOutput;

fileA=argv[1];
fileB= argv[2];
overlapcount =  atoi(argv[3]);
distance= atoi(argv[4]);
fileOutput=argv[5];

if  ((bedA.readTags(fileA))&&(bedB.readTags(fileB)))  {

        bedA.mapTagstoChrom();
        bedB.mapTagstoChrom();

        resultBed=substractBedfiles( &bedA , &bedB,overlapcount,distance);
        resultBed.outputBedTags(fileOutput);
    }else{

          cout << "Error reading data\n";
          cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

    }

    return 0;
}
