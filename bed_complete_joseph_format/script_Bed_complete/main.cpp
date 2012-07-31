#include <iostream>
#include <string>
#include <fstream>
#include  <vector>

using namespace std;


struct BedRead{
 string chrom;
 int chromStart;
 int chromEnd;
 char chromStrand;
};



int main()
{

    string str_filename;
    vector<BedRead> bedfile;
    ifstream i_bedfile;
    BedRead tempRead;
    bool isTag;


    cout<<"Please enter file name from Joseph group for Bed conversation\n";
    cout<<"Output will be [filename]_output.bed\n";

    cin >> str_filename;

    i_bedfile.open(str_filename.c_str() );

    if (!i_bedfile.is_open()){
          cout << "Program failed to open file. Press ENTER to exit";
          cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
          return 0;
    }

    cout<<"1 for tag, 0 for bed\n";
    cin >>isTag;

    while (!i_bedfile.eof())
    {

        i_bedfile >> tempRead.chrom;
        i_bedfile >> tempRead.chromStart;

        if (isTag){
             i_bedfile >> tempRead.chromStrand;
            //Si positif, la fin du read est 36 bp plus loins
            if ( tempRead.chromStrand=='+' )
                tempRead.chromEnd=(tempRead.chromStart+35);
            //Sinon, ChromStart représente en fais la fin du read, donc inverser
            else
                {
                    tempRead.chromEnd=tempRead.chromStart;
                    tempRead.chromStart= (tempRead.chromEnd-35);
                }
        }
        else
        {
            int junk;
             i_bedfile >> junk;

             tempRead.chromEnd=(tempRead.chromStart+25);
             tempRead.chromStart=(tempRead.chromStart-25);
        }


        bedfile.push_back(tempRead);

    }

//Ouput Data
   ofstream o_bedfile;
    str_filename= str_filename+".out";
   o_bedfile.open((str_filename.c_str()));


    for(int i=0; i<(bedfile.size());i++)
    {

     o_bedfile <<(bedfile[i].chrom +"\t");
     o_bedfile <<(bedfile[i].chromStart)<<"\t";
     o_bedfile <<(bedfile[i].chromEnd)<<"\t";

    if (isTag)
        o_bedfile <<(bedfile[i].chromStrand)<<"\n";
     else
        o_bedfile << "\n";
    }


          cout << "Program stop. Press ENTER to exit";
            cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );


return 0;
}


