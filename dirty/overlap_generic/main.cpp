#include <iostream>

using namespace std;


struct fileLine{
std::string chr;
int start;
int end;
string rest;
};


int main(int argc, char* argv[])
{

    string firstArg, secondArg;
    ifstream firstFile, secondFile;
    int overlapNeeded;
    if (argc <=3){
        cerr << "Arguments needed. -h for help";
        return 0;
    }
    firstArg=argv[1];
    secondArg=argv[2];
    overlapNeeded=argv[3];

    firstFile.open(firstArg.c_str());
    secondFile.open(secondArg.c_str());

    vector<fileLine> firstFileVector;
    vector<fileLine> secondFileVector;
    vector<fileLine> weKeep;
    fileLine tempValue;
    char tempstuff[3000],
    if (firstFile.good())&&(secondFile.good()){
        while(!(firstFile.eof()){

                tempValue.chr << firstFile;
                tempValue.start << firstFile;
                tempValue.end <<firstFile;
                tempValue.getline(tempstuff,3000);
                tempValue.rest= tempstuff;
                firstFileVector.push_back(tempValue);
              }

        while(!(secondFile.eof()){

                tempValue.chr << secondFile;
                tempValue.start << secondFile;
                tempValue.end <<secondFile;
                tempValue.getline(tempstuff,3000);
                tempValue.rest= tempstuff;
                secondFileVector.push_back(tempValue);
              }

         vector<fileLine>::iterator firstIt;
         vector<fileLine>::iterator secondIt;





        int overlap=0;

        for(firstIt=firstFileVector.begin(); firstIt!=firstFileVector.end();firstIt++){

            for(secondIt=secondFileVector.begin(); secondIt!=secondFileVector.end();secondIt++){

                if (firstIt->chr==secondIt->chr){

                    overlap=overlapCount(firstIt->start,firstIt->end,secondIt->start, secondIt->end);
                    if (overlap>=overlapNeeded){
                        weKeep.push_back(*firstIt);

                    }

                }
            }
        }



    }
    else


    return 0;
}


//Do these two overlap
bool checkOverlap(int X1, int X2, int Y1, int Y2){

    //Start is inside
    if ((X1>=Y1)&&(X1<=Y2))
        return true;

    //Stop is inside
    if ((X2>=Y1)&&(X2<=Y2))
        return true;

    //Y is englobed by X
    if ((X1<=Y1)&&(X2>=Y2))
            return true;

return false;
}

//Return how many bp overlap
int overlapCount(int X1, int X2, int Y1, int Y2){

  int overlap=0;
  int bigX=X1, smallY=Y1;

    if (checkOverlap(X1, X2, Y1, Y2)){

        if (X2>X1)
            bigX =X2;

        if (Y2<Y1)
            smallY=Y2;
        overlap =smallY- bigX;

        return overlap;
    }
    else
        return overlap;
}




