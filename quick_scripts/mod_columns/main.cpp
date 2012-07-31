#include <iostream>
#include<fstream>
#include <string>
#include "utility.h"
#include <map>
#include <string.h>
using namespace std;

struct sixbed{

string chr;
int start;
int end;
string four;
string five;
string six;

};

struct UCSCID{

string gene;
string id;

};

struct fourbed{

string chr;
int start;
int end;
string ID;

};
using namespace utility;

int main(int argc, char* argv[])
{
ifstream readfile;
string curline, outputID, domainPath,ID, regionPath, genePath,concatPath;

        int count=0;
//        string column[6];
   // ofstream columnoutput( outputID+"1.txt");

    for (int i = 1; i < argc; i++) {
            if (strcmp(argv[i], "-ID")==0) {
                ID = (argv[i+1]);
            } else if (strcmp(argv[i],"-o")==0) {
                outputID = argv[i + 1];
            }else if (strcmp(argv[i],"-d")==0) {
                domainPath = (argv[i+1]);
            }else if (strcmp(argv[i],"-r")==0) {
                regionPath = (argv[i+1]);
            }else if (strcmp(argv[i],"-g")==0) {
                genePath = (argv[i+1]);
            }else if (strcmp(argv[i],"-c")==0) {
                concatPath = (argv[i+1]);
            }
    }

     if ((ID.size()==0)||(outputID.size()==0)||(domainPath.size()==0)||(regionPath.size()==0)){
        cerr<<"Program signature is -g <GeneFilePath> -ID <UniqueID> -o <OutputID> -d <Domain Path> -r <sregion path> -c <concatFile Path>";
       return 0;
    }


        readfile.open(regionPath.c_str());
        if (!(readfile)){
            cerr <<"Invalid file name for culling columns, aborting"<<endl;
            return 0;
        }

     ID +="_";
    ofstream columnoutput( outputID+"1.txt");

    string line;
    while(!std::getline(readfile, line, '\n').eof())
    {
         string chr, start, stop;
        {
            bool result=false;
             utility::Tokenizer ourLine(line);

                /**< Read chr, start and stop positions, abort if there was no token equivalent to */
             result=ourLine.NextToken();
             chr = ourLine.GetToken();
             result=ourLine.NextToken();
             start = ourLine.GetToken();
             result=ourLine.NextToken();
             stop= ourLine.GetToken();

             string OutputID =utility::concatStringInt(ID,count);
          //   if ((chr.size()==0)||(start.size()==0)||(stop.size()==0))
           //     result=false;


             if (result){
                columnoutput << chr<< '\t' << start << '\t' << stop  << '\t'  << OutputID << endl;
                count++;
             }
             else
                {
                    cerr << "Aborting, invalide tab line in " << regionPath.c_str() <<endl;
                     return 0;
                }
        }
    }


  /* while(!std::getline(readfile, column[0], '\t').eof())
    {
        std::getline(readfile, column[1], '\t');
        std::getline(readfile, column[2],'\t');
        //Rest of line
         std::getline(readfile, column[3]);
       // std::getline(readfile, column[3],'\t');


  //      std::getline(readfile, column[5]);

//      Validate this...
       // column [3] =utility::concatStringInt(ID,count);
      string OutputID; =utility::concatStringInt(ID,count);

      //   column[4]=(column[4]+"_");

        columnoutput << column[0]<< '\t' << column[1] << '\t' << column[2]  << '\t'  << OutputID << endl;
        count++;
    } */
     columnoutput.close();
     ofstream overlapoutput(outputID+"2.txt");


     ifstream domaineInput(domainPath);
     ifstream ourRegInput(outputID+"1.txt");
     ifstream geneInput(genePath);

      if (!(domaineInput)){
            cerr <<"Invalid file name for domains, aborting"<<endl;
            return 0;
        }
      if (!(geneInput)){
            cerr <<"Invalid file name for gene ID, aborting"<<endl;
            return 0;
        }

   //read our two datasets input structures.
    cerr << "Reading region data" <<endl;
    map<string, vector<fourbed>> ourRegMap;
    fourbed tempFourBed;
    while(!std::getline(ourRegInput, tempFourBed.chr, '\t').eof())
    {
        string temp;
        tempFourBed.chr=clean_WString(tempFourBed.chr);

        vector<fourbed> * pvector;
        std::getline(ourRegInput, temp, '\t');
        temp=clean_WString(temp);

        tempFourBed.start = utility::stringToInt(temp);
        std::getline(ourRegInput, temp,'\t');
        temp=clean_WString(temp);
        tempFourBed.end = utility::stringToInt(temp);
        std::getline(ourRegInput, tempFourBed.ID);
        tempFourBed.ID=clean_WString(tempFourBed.ID);
        pvector=&(ourRegMap[tempFourBed.chr]);
        pvector->push_back(tempFourBed);
    }

        cerr << "Reading Domain data" <<endl;
    map<string, vector<sixbed>> ourDomainMap;
    sixbed tempSixBed;
 while(!std::getline(domaineInput, tempSixBed.chr, '\t').eof())
    {
        string temp;
        tempSixBed.chr=clean_WString(tempSixBed.chr);
        vector<sixbed> * pvector;
        std::getline(domaineInput, temp, '\t');
         temp=clean_WString(temp);
        tempSixBed.start = utility::stringToInt(temp);
        std::getline(domaineInput, temp,'\t');
        temp=clean_WString(temp);
        tempSixBed.end = utility::stringToInt(temp);
        std::getline(domaineInput, tempSixBed.four,'\t');
        std::getline(domaineInput, tempSixBed.five,'\t');
        std::getline(domaineInput, tempSixBed.six);
        tempSixBed.six=clean_WString(tempSixBed.six);
        pvector=&(ourDomainMap[tempSixBed.chr]);
        pvector->push_back(tempSixBed);
    }
  //   int innercount=0;



cerr << "Reading gene Names" <<endl;
vector<UCSCID> geneID;
map<string,string> geneIDMAP;

UCSCID tempID;
//Read UCSC Data and transform our cluster ID into gene ID
    while(!std::getline(geneInput, tempID.gene, '\t').eof())
    {
        tempID.gene=clean_WString(tempID.gene);
      //  string temp;
        std::getline(geneInput, tempID.id);
         tempID.id=clean_WString( tempID.id);
      //   = utility::stringToInt(temp);
        geneID.push_back(tempID);
      //  cerr << "inserting" <<" "<< tempID.gene <<endl;//<< " " << tempID.id << endl;
     //  cerr << "inserting" <<" "<< tempID.id <<endl;
     //   string a;
     //   cin >> a;
    }


    for (auto USCS: geneID)
        geneIDMAP[USCS.gene]=USCS.id;

cerr << "Finished reading our data" <<endl;
//We overlap both Maps. our output is
//DomainID <List of tab delimited regions overlapping them>
    for(auto it=ourDomainMap.begin(); it!=ourDomainMap.end(); it++){
            vector<sixbed> * psixvector;
            vector<fourbed> * pfourvector;
            psixvector = &(it->second);
            pfourvector = &(ourRegMap[it->first]);
            string output;
            for (auto sixit=psixvector->begin(); sixit!=psixvector->end(); sixit++){
            //We compare each gene/domain to every region and write the ones that overlap to our output string
                    output= sixit->four;
                     //innercount=0;
                        for (auto fourit : *pfourvector){
                          if (utility::isOverlap(sixit->start,sixit->end, fourit.start, fourit.end,OverlapType::OVERLAP_PARTIAL)){
                                output += "\t";
                                output +=fourit.ID;
                            }
                        }
                overlapoutput << output << endl;
            }
    }
    overlapoutput.close();

/*
//Replace the ID with gene ID so we have something to combine
  for(auto it=ourDomainMap.begin(); it!=ourDomainMap.end(); it++){
        //For each ID, replace
     vector<sixbed> * psixvector;
     psixvector = &(it->second);
     for (auto sixit=psixvector->begin(); sixit!=psixvector->end(); sixit++){
             for(auto IDit=geneID.begin(); IDit!=geneID.end(); IDit++){
                    if (IDit->id==strID){
                       // cerr<<"replacing" << endl;
                        strID=IDit->gene;
                        break;
                    }
             }

     }
  }
  cerr << "Finished transforming gene ID" <<endl;
*/
//Read out outputed data and at each read, transform the ID into it's gene name
    map<string, string> overlapData;
    ifstream ourOverlapInput(outputID+"2.txt");
    string tempStr;
    while(!std::getline(ourOverlapInput, tempStr).eof())
    {

        {
        tempStr= clean_WString(tempStr);
        string strID, tempRest;
        utility::Tokenizer data(tempStr);
        data.NextToken();
        strID=data.GetToken();
            for(auto IDit=geneID.begin(); IDit!=geneID.end(); IDit++){
                    if (IDit->id==strID){
                       // cerr<<"replacing" << endl;
                        strID=IDit->gene;
                        break;
                    }
             }
        while (data.NextToken())
            {
                   // string temp;
                    tempRest+="\t";
                    tempRest+=data.GetToken();
                  //  tokens.push_back(data.GetToken());
            }

        overlapData[strID]+=tempRest;
        }
    }
     ofstream overlapoutputunique(outputID+"3.txt");
     ofstream overlacount(outputID+"4.txt");
    //Now we combine similar elements and eliminate duplicates
    map<string,string>::iterator mapit;
    vector<string>::iterator stringit,tokensEnd;
    for(mapit=overlapData.begin(); mapit!=overlapData.end(); mapit++)
    {
        string curRest,ID;
        ID = mapit->first;
        curRest=mapit->second;
        overlapoutputunique << ID;
        overlacount << ID;
        //Tokenize and spit out
            {
                using namespace utility;
                utility::Tokenizer data(curRest);
                vector<string> tokens;
                while (data.NextToken())
                {
                        string temp;
                        temp=data.GetToken();
                        tokens.push_back(data.GetToken());
                }
                 sort(tokens.begin(), tokens.end());
                 tokensEnd=unique(tokens.begin(), tokens.end());

                //Remove duplicates
                        int count=0;
                    //We count the number of unique items
                    for(stringit=tokens.begin(); stringit!=tokensEnd; stringit++){
                        count++;
                        overlapoutputunique <<"\t"<< *stringit ;
                    }
                   // if (tokens.size()>0)
                        overlapoutputunique << "\n";

                 overlacount<< "\t" << count << endl;
                 mapit->second= utility::convertInt(count);
            }
        }


        //Finally, add to our concat file

        ifstream  concatIO(concatPath);

         if (!(concatIO)){
            cerr <<"Invalid file name for concat, aborting"<<endl;
            return 0;
        }

//Line one with header
        string templine, writeresult;
        std::getline(concatIO, templine);

          size_t cr_idx;
          cr_idx = templine.find('\r', 0);
          if (cr_idx != string::npos) {
             templine=(templine.substr(0, cr_idx));
          }


        templine+=("\t"+outputID+"\n");

        writeresult+=templine;

        while(!std::getline(concatIO, tempStr).eof())
        {
           //Get the line, ID and rest.
        tempStr=clean_WString(tempStr);
        string strID, overlapLine, tempRest;
        utility::Tokenizer data(tempStr);
        data.NextToken();
        strID=data.GetToken();
        //cerr << "testing ID " << strID << endl;
        while (data.NextToken())
            {
                   // string temp;
                    tempRest+="\t";
                    tempRest+=data.GetToken();
                  //  cerr << "Next token is" << data.GetToken();
                  //  tokens.push_back(data.GetToken());
            }
       // std::getline(ourOverlapInput, tempRest);
        //Find the completition data and concat

       // utility::debug_string(strID);

        //TEMP
  //   auto curID=geneIDMAP[strID];
       // cerr <<curID <<endl;
        overlapLine=overlapData[strID];

      //  cout <<curID << " "<< overlapLine<<endl;
        if (overlapLine.size()==0)
            overlapLine ="0";
        writeresult+=strID;
        writeresult+=(tempRest);
        writeresult+=("\t"+overlapLine+"\n");
             //Our data
           // overlapData
        }
        concatIO.close();
        ofstream  concatWrite(concatPath);
        concatWrite<<writeresult;

}



