#include <bedFile.h>
#include "BedUtility.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>


using namespace std;


int main(){


bedFile inputData, comparisonData;

barrierResults resultData;
string getChromfromint(int chromNum);

string str_filename, str_comparison;
int windowSize, tagperIslands, tagperWindows, nbemptyWindows;

//Read everything from a bedFile.

   if(!inputData.openReadTags())
        {
        cout << "Program failed to open file. Press ENTER to exit";
        cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

        return 0;
        }

        cout << "Please input window size, tag per islands, tag per window and nb of empty windows.\n";
        cin >> windowSize;
        cin >> tagperIslands;
        cin >> tagperWindows;
        cin >> nbemptyWindows;

   if(windowSize<50)
        {
        cout << "Invalide window size. Press ENTER to exit";
        cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        return 0;
        }


    //Map tags to windows.
    inputData.mapTagstoWindows(windowSize);

    //Generate the islands, for now default 8.
    inputData.countIslands(tagperIslands,nbemptyWindows,tagperWindows);


    //Now we want to overlap our data...
    // Read our protein data.
    cout<<"Please enter file name containing data we want to compare with our islands.\n";
    cout<<"File must contain three columns, chrom, chromStart and chromEnd\n";
    cin >> str_comparison;

    //Read the bedFile data
   if (!comparisonData.readTags(str_comparison)){
       return 0;
   }

    comparisonData.mapTagstoChrom();

    //Now we check if protein from comparisonData is near the frontiers of islands from inputData
   resultData= barrierOverlap( inputData, comparisonData, 1000) ;



    ofstream o_Bar, o_Inside, o_Others, o_Islands;


    o_Bar.open("barrier.txt");
    o_Inside.open("inside.txt");
    o_Others.open("other.txt");
    o_Islands.open("islands.txt");

    for (int i=0; i<CHROMCOUNT; i++){
         ChromSum *pCurChrom;
        pCurChrom =  inputData.getChrom(i);
//Barrier results
       for (unsigned int k=0; k <resultData.barrierTags[i].size();k++){
        o_Bar << getChromfromint(resultData.barrierTags[i].at(k).chrom) << "\t";
        o_Bar << resultData.barrierTags[i].at(k).chromStart<< "\t";
        o_Bar << resultData.barrierTags[i].at(k).chromEnd<<"\n";
      }

    for (unsigned int k=0; k <resultData.centerTags[i].size();k++){
//Islands results
        o_Inside << getChromfromint(resultData.centerTags[i].at(k).chrom) << "\t";
        o_Inside << resultData.centerTags[i].at(k).chromStart<< "\t";
        o_Inside << resultData.centerTags[i].at(k).chromEnd<<"\n";
      }
//Other tags
       for (unsigned int k=0; k <resultData.otherTags[i].size();k++){
        //Output each result;
        o_Others << getChromfromint(resultData.otherTags[i].at(k).chrom) << "\t";
        o_Others << resultData.otherTags[i].at(k).chromStart<< "\t";
        o_Others << resultData.otherTags[i].at(k).chromEnd<<"\n";
      }

       for (unsigned int k=0; k <pCurChrom->iList.size();k++){


        //Output each Island;
        if (pCurChrom->iList.at(k).iWindowCount>1){
            o_Islands << getChromfromint(pCurChrom->iList.at(k).iChrom) << "\t";
            o_Islands << pCurChrom->iList.at(k).iStart<< "\t";
            o_Islands << pCurChrom->iList.at(k).iStop<< "\t";
            o_Islands << pCurChrom->iList.at(k).iWindowCount << "\t";
            o_Islands << pCurChrom->iList.at(k).iCount<<"\n";

        }
        }



    }

  cout << "Program finished. See output text files for results. Press ENTER to exit";
        cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

return 0;
}


