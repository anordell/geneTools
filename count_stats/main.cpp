#include "bedFile.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <dirent.h>

using namespace std;

string GetFileExtension(const std::string& FileName);
vector<string> listFile();
std::string convertInt(int number);

int main(int argc, char* argv[])
{
    vector<string> list_files;
    int tempcount;


    if (argc!=2)
        return 0;

    list_files=listFile();
    bedFile tempBed;
    string results;
    ofstream o_output;

    for (unsigned int i=0;  i<list_files.size() ; i++  ){

        tempBed.reset();
        tempBed.readTags(list_files.at(i));
        results += list_files.at(i);
        results += " : ";
        tempcount=tempBed.tagCount();
       // if (tempcount>1)
      //      tempcount--;
        results += convertInt(tempcount) ;
        results +="\n";
    }

    o_output.open(argv[1]);
    o_output << results;

    return 0;
}


vector<string> listFile(){
        vector<string> results;
        DIR *pDIR;
        struct dirent *entry;
        if( pDIR=opendir(".") ){
                while(entry = readdir(pDIR)){
                        if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 ){

                       // cout<<(GetFileExtension(entry->d_name));
                            if ((GetFileExtension(entry->d_name)=="bed")||(GetFileExtension(entry->d_name)=="txt"))
                                 results.push_back(entry->d_name);
                        }

                }
                closedir(pDIR);
        }

    return results;
}


string GetFileExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}


//Utility function courtesy of the internet.
string convertInt(int number)
{
    if (number == 0)
        return "0";
    string temp="";
    string returnvalue="";
    while (number>0)
    {
        temp+=number%10+48;
        number/=10;
    }
    for (unsigned int i=0;i<temp.length();i++)
        returnvalue+=temp[temp.length()-i-1];
    return returnvalue;
}

