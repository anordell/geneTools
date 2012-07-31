//######################################
//Program by Alexei Nordell Markovits (2011)
//Alexei.Nordell@usherbrooke.ca
//
//All rights reserved. You may not distribute or modify this program for purpose of distribution without the written
//consent of it's creator
//#######################################



#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdlib.h>

using namespace std;


float totalcount=0;

struct nucleo_data{
string chrom;
int pos;
char base;
long int count_sens;
long int count_anti_sens;
long int snp_maj;
long int snp_min;


};


struct gene{
string id;
string locus_tag;
string alt_locus_tag;
string    gi;
char strand;
long int begin;
long int end;
string name;
string desc;
string cog;
string ec;
string num;
string evidence;
string restofit;



int count_sens;
int count_anti_sens;
int snp_maj;
int snp_min;


};


vector<nucleo_data> readPileup(ifstream& input);

vector<gene> readGenbank(ifstream& input);
int findLastEnd(vector<gene> allgenes);
vector<gene> generateIntergenic(vector<gene>);
string concatStringInt(string ourstring, int ourInt);

int main(int argc, char* argv[])

{

    float ONEBILLION=1000000000;
    string tagFileName;
    string geneFileName;
    string outputFileName,outputTagName;
    ifstream i_bedfile;
    vector<nucleo_data> allreads;
    vector<gene> allgenes;
    nucleo_data tempread;


    if (argc < 5){
        cout<< " Signature is  <tagfile> <genomefile> <Geneoutputname> <PileupReportOutput>";
        return 0;
    }
    tagFileName=argv[1];
    geneFileName=argv[2];
    outputFileName=argv[3];
    outputTagName=argv[4];


    i_bedfile.open(tagFileName.c_str() );
    if (!i_bedfile.is_open()){
          return 0;
    }


 allreads= readPileup(i_bedfile);

i_bedfile.close();

//Read Gene data
i_bedfile.clear();

i_bedfile.open(geneFileName.c_str() );
    if (!i_bedfile.is_open()){
          return 0;
    }

    allgenes=readGenbank(i_bedfile);
  //  int gene_size= findLastEnd(allgenes);

    i_bedfile.close();



    allgenes=generateIntergenic(allgenes);


//Count tags per gene
//For gene
 for (unsigned int k=0; k <allgenes.size();k++){
    //For each nuclotide in the gene
    //If gene inverse, we reverse.
    int start,end, buffer;
    start= allgenes.at(k).begin;
    end=allgenes.at(k).end;
    if (end<start){
        buffer=end;
        end=start;
        start=buffer;
    }

     for (int i=(start-1); i <= (end-1); i++){
         allgenes.at(k).count_sens+=allreads.at(i).count_sens;
         allgenes.at(k).count_anti_sens+=allreads.at(i).count_anti_sens;
         allgenes.at(k).snp_maj+=allreads.at(i).snp_maj;
         allgenes.at(k).snp_min+=allreads.at(i).snp_min;
    }

 }

    float mult;
   mult=ONEBILLION/totalcount;
    int buffer;


    ofstream o_output;
    o_output.open((outputFileName.c_str()));
    int size;
    float count;
//For each gene
       for (unsigned int k=0; k <allgenes.size();k++){
             //Calc our stats


                if ( allgenes.at(k).strand=='-'){
                    buffer=allgenes.at(k).begin;
                    allgenes.at(k).begin=allgenes.at(k).end;
                    allgenes.at(k).end=buffer;
          }

                o_output << allgenes.at(k).id << "\t";
                 o_output << allgenes.at(k).locus_tag << "\t";
                 o_output << allgenes.at(k).alt_locus_tag << "\t";
                o_output << allgenes.at(k).name << "\t";
                 o_output << allgenes.at(k).desc << "\t";

                 if (allgenes.at(k).strand=='-')
                    size=(allgenes.at(k).begin- allgenes.at(k).end);
                else
                  size=(allgenes.at(k).end-allgenes.at(k).begin);

                size++;
                o_output << size<< "\t";
                o_output << allgenes.at(k).begin << "\t";
                o_output << allgenes.at(k).end << "\t";
                 o_output << allgenes.at(k).strand << "\t";

                o_output << allgenes.at(k).count_sens << "\t";
//  Divide and normalise
                count= allgenes.at(k).count_sens;
                if (count!=0){
                    count= count/size;
                    count=count*mult;
                }
                o_output << count << "\t";

                 o_output << allgenes.at(k).count_anti_sens << "\t";

                count= allgenes.at(k).count_anti_sens;
               if (count!=0){
                    count= count/size;
                    count=count*mult;
               }
                o_output << count << "\n";

            }


    o_output.close();
    o_output.open((outputTagName.c_str()));


       for (unsigned int k=0; k <allreads.size();k++){
                o_output << allreads.at(k).chrom << "\t";
                o_output << allreads.at(k).pos << "\t";
                o_output << allreads.at(k).base << "\t";
                o_output << allreads.at(k).count_sens << "\t";
                o_output << allreads.at(k).count_anti_sens << "\t";
                o_output << allreads.at(k).snp_maj << "\t";
                o_output << allreads.at(k).snp_min << "\n";

            }

    return 0;
}







vector<nucleo_data> readPileup(ifstream& input){


    vector<nucleo_data> allreads;
    nucleo_data tempread;
    int depth;
    char buffer;

while (!input.eof())
    {
        input >> tempread.chrom;
         if( input.eof() ) break;

            tempread.count_sens=0;
            tempread.count_anti_sens=0;
            tempread.snp_maj=0;
            tempread.snp_min=0;


           input>>tempread.pos;
           input>>tempread.base;
           input>> depth;

            while (buffer!='@')
                buffer=input.get();


         while (buffer!='\n'){
                if( input.eof() ) break;
                buffer=input.get();

             switch(buffer)
             {
              case ',':
                 tempread.count_sens++;
                 totalcount++;
                 break;
              case '.':
                tempread.count_anti_sens++;
                totalcount++;
              break;
                }

            if ((buffer>='a')&&(buffer<='z')){
                tempread.snp_min++;
                totalcount++;
            }

            if ((buffer>='A')&&(buffer<='Z')){
                tempread.snp_maj++;

            }
         }

         allreads.push_back(tempread);
            //Ignore rest of line
           // i_bedfile.ignore(5000, '\n');
    }
    input.close();

return allreads;

}


 vector<gene>readGenbank(ifstream& input){



 char tempdata[2500];
  vector<gene> allgenes;
  gene tempgene;
  int buffer;

    //Discard header
    input.ignore(5000, '\n');

  while (!input.eof())
    {

        input.getline(tempdata , 2500, '\t');
         if( input.eof() )
             break;
        tempgene.id.assign(tempdata);

        input.getline(tempdata , 2500, '\t');
        tempgene.locus_tag.assign(tempdata);

        input.getline(tempdata , 2500, '\t');
        tempgene.alt_locus_tag.assign(tempdata);

        input.getline(tempdata , 2500, '\t');
        //tempgene.gi=atoi(tempdata);
        tempgene.gi.assign(tempdata);

        input.getline(tempdata , 2500, '\t');
        tempgene.strand=(tempdata[0]);

        input.getline(tempdata , 2500, '\t');
        tempgene.begin= atoi(tempdata);

         input.getline(tempdata , 2500, '\t');
        tempgene.end= atoi(tempdata);

        input.getline(tempdata , 2500, '\t');
        tempgene.name.assign(tempdata);

         input.getline(tempdata , 2500, '\t');
         tempgene.id.assign(tempdata);

//Store the rest the rest
        input.getline(tempdata,2500, '\n');
        tempgene.restofit.assign(tempdata);

          if (tempgene.strand=='-'){
                buffer=tempgene.begin;
                tempgene.begin=tempgene.end;
                tempgene.end=buffer;
          }


        tempgene.count_anti_sens=0;
        tempgene.count_sens=0;
        tempgene.snp_maj=0;
        tempgene.snp_min=0;

        allgenes.push_back(tempgene);

        if( input.eof() )
         break;

    }

return allgenes;
 }


int  findLastEnd( vector<gene> allgenes){

int last=0;

    for (unsigned int k=0; k <allgenes.size();k++){
        //Find last gene
        if (allgenes.at(k).end>last)
            last=allgenes.at(k).end;
    }

return last;
}

vector<gene> generateIntergenic(vector<gene> allgenes){

gene currentgene, lastgene, intergenic;
int regionCount=0;
vector<gene> completedGenes;
string regionName;

if(allgenes.at(0).begin!=1){
    regionName= concatStringInt("IntergenicRegion_",regionCount);
    regionCount++;
    intergenic.id=regionName;
    intergenic.begin=1;
    intergenic.end=allgenes.at(0).begin-1;
    intergenic.strand='+';

     intergenic.count_anti_sens=0;
    intergenic.count_sens=0;
    intergenic.snp_maj=0;
    intergenic.snp_min=0;
    completedGenes.push_back(intergenic);
}
lastgene=allgenes.at(0);

completedGenes.push_back(lastgene);


    for (unsigned int k=1; k <(allgenes.size());k++){

        currentgene=allgenes.at(k);

        //If last gene ends before current gene begins, add an intergenic region
        if ((lastgene.end+1)<currentgene.begin)
        {
            regionName= ("IntergenicRegion_"+lastgene.locus_tag);
            if (k<(allgenes.size()-1) )
                regionName=regionName+allgenes.at(k).locus_tag;
            regionCount++;
            intergenic.id=regionName;
            intergenic.begin=lastgene.end+1;
            intergenic.end=currentgene.begin-1;
            intergenic.strand='+';

             intergenic.count_anti_sens=0;
            intergenic.count_sens=0;
            intergenic.snp_maj=0;
            intergenic.snp_min=0;

            completedGenes.push_back(intergenic);
        }
        //add current gene
            completedGenes.push_back(currentgene);
        //swap
            lastgene=currentgene;
    }
return completedGenes;
}


string concatStringInt(string ourstring, int ourInt){

    std::stringstream sstm;
    sstm << ourstring << ourInt;
    return  sstm.str();

}

