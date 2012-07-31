
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include "main.h"

using namespace std;

int main(int argc, char* argv[])
{
     ifstream inputStream;
     ofstream outputStream;

    string inputpath, outputpath;
        for (int i = 1; i < argc; i++) {
                if (strcmp(argv[i], "-i")==0) {
                    // We know the next argument *should* be the filename:
                    inputpath = argv[i + 1];
                } else if (strcmp(argv[i],"-o")==0) {
                    outputpath = argv[i + 1];
                }
        }

 if ((inputpath.size()==0)||(outputpath.size()==0)){
        cerr<<"Program signature is -i <inputpath> -o <outputpath>";
        return 0;
    }

  inputStream.open(inputpath.c_str());
        if (inputStream.good()){
        }
        else{
            cerr << "Error in opening inputpath";
            return 0;
        }
      outputStream.open(outputpath.c_str());
        if (outputStream.good()){
        }
        else{
            cerr << "Error in opening outputpath";
            return 0;
        }

        map<string,string> ourMap;

//Iterator over our file.
while(!(inputStream.eof())){
        {

        string ID, rest;
        inputStream >> ID;
        getline(inputStream,rest);
        ourMap[ID]+=rest;

        cerr << "I am inserting ID " << ID << endl;
       // cerr << " And creating Key " << ourMap[ID];
        inputStream.peek();
        if (inputStream.eof())
            break;
        }

      }

    map<string,string>::iterator mapit;
    vector<string>::iterator stringit,tokensEnd;
    for(mapit=ourMap.begin(); mapit!=ourMap.end(); mapit++)
    {
        string curRest,ID;
        ID = mapit->first;
        curRest=mapit->second;
        outputStream << ID <<"\t";
        //Tokenize and spit out
            {
                Tokenizer data(curRest);
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
                    for(stringit=tokens.begin(); stringit!=tokensEnd; stringit++){

                        outputStream << *stringit << "\t";
                    }
                    if (tokens.size()>0)
                        outputStream << "\n";
            }
    }
}

//CPP file
//Courtesy of the internet, thank you mister anonymous at
//http://stackoverflow.com/questions/53849/how-do-i-tokenize-a-string-in-c
const string Tokenizer::DELIMITERS(" \t\n\r");

Tokenizer::Tokenizer(const std::string& s) :
    m_string(s),
    m_offset(0),
    m_delimiters(DELIMITERS) {}

Tokenizer::Tokenizer(const std::string& s, const std::string& delimiters) :
    m_string(s),
    m_offset(0),
    m_delimiters(delimiters) {}

bool Tokenizer::NextToken()
{
    return NextToken(m_delimiters);
}

bool Tokenizer::NextToken(const std::string& delimiters)
{
    size_t i = m_string.find_first_not_of(delimiters, m_offset);
    if (string::npos == i)
    {
        m_offset = m_string.length();
        return false;
    }

    size_t j = m_string.find_first_of(delimiters, i);
    if (string::npos == j)
    {
        m_token = m_string.substr(i);
        m_offset = m_string.length();
        return true;
    }

    m_token = m_string.substr(i, j - i);
    m_offset = j;
    return true;
}
 const std::string Tokenizer::GetToken() const
{

    return m_token;

}
