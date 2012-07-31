#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>


using namespace std;




struct SNP{
string name;
int pos;
char base;
char basemodif;
double numbers[5];
char letter1;
int number2;
char letter2;
char restOfLine[600];

};

typedef map<int, SNP> mapofSNP;

void outputSnpfile(ofstream& output, mapofSNP& outputMap);
void readSnpFile(ifstream& input, mapofSNP& inputMap);
mapofSNP findUniquefromA(mapofSNP& inputMap, mapofSNP& compareMap);
mapofSNP findCommon(mapofSNP& inputMap, mapofSNP& inputMap2);


int main(int argc, char* argv[])
{
    string file1Name, file2Name, output1, output2, outputcommon;

    ifstream isFile1, isFile2;
    ofstream osFile1, osFile2, osFileCommon;

     mapofSNP mapSnpFile1;
     mapofSNP mapSnpFile2;

     mapofSNP mapSnpUnique1;
     mapofSNP mapSnpUnique2;
     mapofSNP mapSnpCommon;

    if (argc < 6){
        cout<< " Signature is  <input 1> <input 2> <outputunique1> <outputunique2> <commonoutput>";
        return 0;
    }
    file1Name=argv[1];
    file2Name=argv[2];
    output1=argv[3];
    output2=argv[4];
    outputcommon=argv[5];

    osFile1.open(output1.c_str());
    osFile2.open(output2.c_str());
    osFileCommon.open(outputcommon.c_str());

    isFile1.open(file1Name.c_str());

    isFile2.open(file2Name.c_str());

    readSnpFile(isFile1, mapSnpFile1);
    readSnpFile(isFile2, mapSnpFile2);


    mapSnpUnique1=findUniquefromA(mapSnpFile1,mapSnpFile2);

    mapSnpUnique2=findUniquefromA(mapSnpFile2,mapSnpFile1);

    mapSnpCommon=findCommon(mapSnpFile1,mapSnpFile2);

    outputSnpfile(osFile1, mapSnpUnique1);
    outputSnpfile(osFile2, mapSnpUnique2);
    outputSnpfile(osFileCommon, mapSnpCommon);

    return 0;
}

mapofSNP findUniquefromA(mapofSNP& inputMap, mapofSNP& compareMap){


      mapofSNP returnMap;
      mapofSNP::iterator iterMap;
      SNP tempSNP;

      int count;
      bool unique=false;

      count =inputMap.size();

     for (iterMap = inputMap.begin(); iterMap != inputMap.end(); ++iterMap) {
        unique=false;

      if  (compareMap.count(iterMap->first)==0)
        {
            unique=true;
        }else
          {
              tempSNP = compareMap[iterMap->first];

              if (tempSNP.basemodif!=iterMap->second.basemodif)
                  unique=true;
          }


     if (unique)
         returnMap.insert(pair<int, SNP>(iterMap->first, iterMap->second));
     }

return returnMap;
}

mapofSNP findCommon(mapofSNP& inputMap, mapofSNP& inputMap2){

    SNP tempSNP;
      mapofSNP returnMap;
      mapofSNP::iterator iterMap;

      bool common=false;

     for (iterMap = inputMap.begin(); iterMap != inputMap.end(); ++iterMap) {

       common=false;

      if  (inputMap2.count(iterMap->first)==1){

         tempSNP = inputMap2[iterMap->first];

              if (tempSNP.basemodif==iterMap->second.basemodif)
                  common=true;

      }

      if (common)
        returnMap.insert(pair<int, SNP>(iterMap->first, iterMap->second));

     }

return returnMap;



}




void readSnpFile(ifstream& input, mapofSNP& inputMap){

    SNP tempSnp;
    char temp;

    temp= input.peek();

    if (temp=='>')
        input.ignore(1000, '\n');
    else
  while (!input.eof())
    {
      //  input >> tempSnp.name;
        input >> tempSnp.pos;
         if( input.eof() ) break;


        input >> tempSnp.base;
        input >> tempSnp.basemodif;


        input.getline(tempSnp.restOfLine, 600);

    /*    input >> tempSnp.numbers[0];
        input >> tempSnp.numbers[1];
        input >> tempSnp.numbers[2];
        input >> tempSnp.numbers[3];
        input >> tempSnp.numbers[4];


        input >> tempSnp.letter1;
        input >> tempSnp.number2;
        input >> tempSnp.letter2;*/

        inputMap.insert(pair<int, SNP>(tempSnp.pos, tempSnp));

    }

    return;
}

void outputSnpfile(ofstream& output, mapofSNP& outputMap){

    SNP tempSNP;

 //Typedef used
    mapofSNP::iterator iterMap;

     for (iterMap = outputMap.begin(); iterMap != outputMap.end(); ++iterMap) {
        tempSNP=(iterMap->second);

      //  output << tempSNP.name << "\t";
        output <<tempSNP.pos << "\t";
        output <<tempSNP.base << "\t";
        output << tempSNP.basemodif<< "\t";
        output << tempSNP.restOfLine << "\n";


       /* output << tempSNP.numbers[0]<< "\t";
        output << tempSNP.numbers[1]<< "\t";
        output << tempSNP.numbers[2]<< "\t";
        output << tempSNP.numbers[3]<< "\t";
        output << tempSNP.numbers[4]<< "\t";


        output << tempSNP.letter1<< "\t";
        output << tempSNP.number2<< "\t";
        output << tempSNP.letter2<< "\n";*/
    }
    return;
}

