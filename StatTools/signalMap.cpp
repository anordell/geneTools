#include "signalMap.h"
#include <tclap/CmdLine.h>
#include "utility.h"
#include <map>
#include "gnuplot_i.hpp"
using namespace std;


class Gnuplot;

void signalMap(int argc, char* argv[])
{
    try
    {
        TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
        /**< Declare and add arguments */
        TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : peakStat", false," ","");
        TCLAP::MultiArg<std::string> fileArg("f","file","List of files containing genomic signals",true,"path to file",cmd);
        TCLAP::MultiArg<std::string> nameArg("n","name","Name of each condition associated with a -f parameter.",true,"ID",cmd);
        TCLAP::ValueArg<std::string> outputArg("o","output","Path to output file",false,"","path to output file",cmd);

        cmd.add(command);
        /**< Parse */
        cmd.parse( argc, argv );



        /**< Assign */
        vector<std::string> filePathVector = fileArg.getValue();
        vector<std::string> nameVector= nameArg.getValue();
        std::string outputPath = outputArg.getValue();
       // ifstream fileStream;
      //  utility::loadStream(filePath, fileStream);


        if (filePathVector.size()!=nameVector.size())
            throw 10;



         /**< Main data structure */
           map<std::string, map<std::string, std::vector<float>>> signalMap;

        int i=0;
        for (string ID: nameVector)
        {
            {
            ifstream inputStream;
            utility::loadStream(filePathVector.at(i),inputStream);
            signalMap[ID]= loadSignalDataFromFile(inputStream);
            i++;
            }
        }


/**< All our data is loaded */
/**< For each map, take every signal associated with the ID in the FIRST signal map and create a heatmap with it. */
/**<  Again, we assume all ID are present for all */

        auto refMap=signalMap.begin();
        for (auto IdMap:refMap->second)
        {
            auto ID = IdMap.first;
            vector<vector<float>> signalVector;
            for (auto currentSignalMap: signalMap )
            {
                if (currentSignalMap.second.count(ID))
                    signalVector.push_back(currentSignalMap.second[ID]);
                else
                {
                    cerr << "missing data, ID " <<ID <<endl;
                    abort();
                }
            }
           auto distData =  generateDistMatrice(signalVector);
           writeHeatMaptoPS(distData, outputPath+ID,nameVector,ID);
        }

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (int e)
    {
        if (e==10)
            std::cerr << "There must be an equal amount of -n and -f parameters. "<< std::endl;
    }
}

/** \brief Write our distance matrix to a GnuPlot ps Heatmap
 *
 * \param heatMap vector<vector<float>> distance Matrix, must be N by N, not verified by program
 * \param psName const string& Filoename to output to
 * \param const string & titlename = "null" Title of graph
 * \return void
 *
 */
void writeHeatMaptoPS(const vector<vector<float>> & heatMap, const string & psName,const vector<string> & legendVector, const string & titlename )
{

    /**< Write directly to PS */

    Gnuplot g1("Heatmap");
    g1.savetops(psName);
    float maxElem=0;
    string gnuCommand;
    float range = heatMap.size()-0.5;
    float xPos=0;
    float yPos=0;

    for(auto column:heatMap)
        maxElem = max( *max_element(column.begin(),column.end()),  maxElem);

    gnuCommand+= "set title \""+titlename+"\"\n";
    gnuCommand+= "unset key\n";
    gnuCommand+= "set tic scale 0\n";
    gnuCommand+= "set palette rgbformula 21,22,23\n";
    gnuCommand+= "set cbrange [0:"+utility::numberToString(maxElem)+"]\n";
    gnuCommand+= "set cblabel \"Score\"\n";
    gnuCommand+= "unset cbtics\n";
    /**< Use this as an example of our legend */
    for (auto legend: legendVector){
        gnuCommand+= "set label \""+legend+"\" at "+utility::numberToString(xPos)+",-1\n";
        xPos+=1;
        }

    for (auto legend: legendVector){
        gnuCommand+= "set label \""+legend+"\" at -1,"+utility::numberToString(yPos)+"\n";
        yPos+=1;
        }

    gnuCommand+= "set xrange [-0.5:"+utility::numberToString(range)+"]\n";
    gnuCommand+= "set yrange [-0.5:"+utility::numberToString(range)+"]\n";
    gnuCommand+=  "set view map\n";
    gnuCommand+="splot '-' matrix with image\n";

    for (auto column: heatMap)
    {
         for (auto row: column)
         {
            gnuCommand+=utility::numberToString(row)+" ";
            cerr<< utility::numberToString(row)+" ";
        }
        gnuCommand+="\n";
        cerr<<"\n";
    }
    gnuCommand+="e\n";
    gnuCommand+="e\n";

    utility::pause_input();


    g1.cmd(gnuCommand);


}


/** \brief Load any number of signals from a single file.
 *
 * \return map<string,vector<float>>  Map with string being name of a signal Type and vector<float> being the signal
 *
 */
std::map<string,vector<float>> loadSignalDataFromFile(ifstream & readfile){

    string line;
   // ifstream readline;
    string curID;
    map<string,vector<float>> returnMap;
    while(!std::getline(readfile, line, '\n').eof())
    {
        /**< Read one line on two for title, we assume file is correctly formated */
        if(curID.size()==0)
            curID=line;
        else
        {
                using namespace utility;
                Tokenizer data(line);
                vector<float> tokens;
                while (data.NextToken())
                {
                        string temp;
                        temp=data.GetToken();
                        tokens.push_back(stringToNumber(temp));
                }
                returnMap[curID]=tokens;
                curID.clear();
        }

    }
    cerr << "Loaded "<< returnMap.size() << " different Signals" <<endl;
return returnMap;
}


vector<vector<float>> generateDistMatrice(vector<vector<float>> listOfSignals){


vector<vector<float>> returnMatrix;




returnMatrix.resize(listOfSignals.size());
/**< Correct size */
for (auto & elem:returnMatrix )
    elem.resize(listOfSignals.size());

    for(int x=0; x< (int)listOfSignals.size(); x++)
    {
        for(int y=0; y< (int)listOfSignals.size(); y++){
          //      cerr <<"Access to X and Y " << x << " " << y <<endl;

           /*     cout << "X is" <<endl;
         for(auto dataX:listOfSignals.at(x) )       {
                cout << dataX << " ";
                }
                cout <<endl;

         cout << "Y is" <<endl;
         for(auto dataY:listOfSignals.at(y) )       {
                cout << dataY << " ";
                }
                cout <<endl; */

         returnMatrix[x][y]= clustering::align_distance(listOfSignals.at(x),listOfSignals.at(y),0,listOfSignals.at(x).size()/10);
        // cerr <<"Result is " <<returnMatrix[x][y]<<endl;
        }
    }

return returnMatrix;
}
