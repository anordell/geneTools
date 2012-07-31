#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>


using namespace std;





int main(int argc, char* argv[])
{

    string inputname, outputname;
    ofstream outputfile;
    ifstream inputfile;
    string chrom;
    int count;

    // outputBed.open("all.bed");

    inputfile.open("C:\\Work\\data\\ours\\H2AZ\\H2AZq10.sam.bed.stats");

    int counting[11]={0,0,0,0,0,0,0,0,0,0,0 };
       while(!(inputfile.eof())){

        inputfile >>chrom ;
         inputfile>> count  ;

        if (count>=10)
            counting[10]++;
        else
            counting[count]++;
    }

    for (int i=0; i<11;i++)
        cout <<counting[i] << "\n";



    return 0;
}
