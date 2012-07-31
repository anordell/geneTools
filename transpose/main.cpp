#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

int main(int argc, char* argv[])
{

    string input, output;

   // input = argv[1];
   // output = argv[2];

    cout << "Nom du fichier"<<"\n";
    cin >> input;
    cout <<"Nom du fichier sortie"<<"\n";
    cin>> output;



  ifstream i_data;
  vector<char>  data;

    i_data.open(input.c_str() );
    if (!i_data.is_open()){
          return false;
    }

 char tempdata[250];
 char count;

  i_data.getline(tempdata , 150, '\t');
    i_data.getline(tempdata , 150, '\t');
      i_data.getline(tempdata , 150, '\n');
  while (!i_data.eof())
    {
          i_data.getline(tempdata , 150, '\t');
         if( i_data.eof() ) break;

       i_data.getline(tempdata , 150, '\t');
       //i_data.getline(count , 1, '\n');
       i_data >> count;
       data.push_back(count);
    }
    i_data.close();


     ofstream o_output;
    o_output.open(output.c_str());



       for (unsigned int k=0; k <data.size();k++){
                o_output << data.at(k) << " ";

            }


    o_output.close();


    }
