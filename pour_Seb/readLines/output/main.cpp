#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>


using namespace std;





int main(int argc, char* argv[])
{

    string inputname, outputname, currentline,currentQueryline, columnStr, lineStr;
    int curValue=0, totalEmpty=0;
    ofstream outputfile;
    ifstream inputfile;
    char temp_char[1000];

    inputname=argv[1];

    inputfile.open(inputname.c_str());

    inputfile.getline(temp_char,1000);
    currentQueryline=temp_char;
       currentline=currentQueryline;
    while (inputfile.eof()==false){
        inputfile.getline(temp_char,1000);
        currentline=temp_char;
        if ( currentline.find("Query:")!=-1){
            if (curValue==0)
                totalEmpty++;
         //   cout <<currentQueryline<< " " << curValue <<"\n";
            curValue=0;
            currentQueryline= currentline;
        }
        else
            curValue++;
    }
    if (curValue==0){
        //cout <<currentQueryline<< " " << curValue <<"\n";
        totalEmpty++;
        //currentQueryline= currentline;

    }
    cout << totalEmpty;

    return 0;
}
