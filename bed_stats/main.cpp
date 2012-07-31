#include "bedFile.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

int main(int argc, char* argv[])
{
 bedFile bedA,bedB, resultBed;
 ofstream filewrite;

int centralPos,curDist;

vector<int> list_centralpos[25];
vector<int> distances;

string fileA, fileB, fileOutput;

fileA=argv[1];
fileOutput=argv[2];
    if ( bedA.readTags(fileA)){
        bedA.mapTagstoChrom();
        bedA.sortChromTags();

        //Found central position of all.
        for (int i=0; i<CHROMCOUNT; i++)
             {
                    ChromSum *pChromA;
                    pChromA = bedA.getChrom(i);

                   for (unsigned int k=0; k <pChromA->tList.size();k++)
                   {

                        centralPos= ((pChromA->tList.at(k).chromStart+pChromA->tList.at(k).chromEnd)/2);
                        list_centralpos[i].push_back(centralPos);
                   }


                  for(unsigned int z=0; z< (list_centralpos[i].size()-1);z++){

                        curDist=list_centralpos[i].at(z+1)-list_centralpos[i].at(z);

                      distances.push_back(curDist);
                  }

             }


        filewrite.open(fileOutput.c_str());

        for(unsigned n=0; n<distances.size();n++)
            filewrite<<distances.at(n);



    }else{

          cout << "Error reading data\n";
          cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
    }


    return 0;
}

