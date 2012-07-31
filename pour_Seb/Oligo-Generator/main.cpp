#include <iostream>
#include <string>
#include <vector>
using namespace std;

struct pGraphs{

 vector<int> nodes;

};


struct oligoNode{
    string composition;

  //  vector<pDist> ED;
  //vector<*oligoNode> pNodes;
};


vector<string> generateNucleoSubgroups(int index, vector<string> nucleo);
int editDist( string curString, string otherString);
bool isClique(pGraphs *ourClique, const vector<oligoNode> &oligoNodes, int newNode , const int dist);
vector<pGraphs> generateNucleoClique(int index, const vector<oligoNode> &oligoNodes, vector<pGraphs> oligoCliques , const int dist );
int editCount(const pGraphs &ourClique,const vector<oligoNode> &oligoNodes);

pGraphs generateNucleoCliqueBetter(int index, const vector<oligoNode> oligoNodes, pGraphs oligoCliques , const int dist );

int totalCount=0;

int main()
{
    int sgSize, gSize;

    //const int NBNUCLEO=6;
   // const int MINEDIT =2;
    const int nucCount=16;
    const int minEditDist=2;

    vector<string> oligoStrings;
    sgSize=3;
  //  gSize=9;

    oligoStrings=generateNucleoSubgroups(sgSize,oligoStrings);

cerr << "There are" << oligoStrings.size() << " strings" <<endl;
    //Generate every possibility of our group.


    //Create our basic graph structure
    vector<oligoNode> oligoGraph;
     for (unsigned int k=0; k< oligoStrings.size();k++){
       oligoNode tempNode;
       tempNode.composition= oligoStrings.at(k);
       oligoGraph.push_back(tempNode);
    }


pGraphs tempGraph;

    tempGraph=  generateNucleoCliqueBetter(nucCount, oligoGraph, tempGraph, minEditDist );

    for( int k =0 ;k< tempGraph.nodes.size(); k++ ){
        cout <<  oligoGraph.at(tempGraph.nodes.at(k)).composition  << " " ;

    }




/*
//Compute all distances
//Could use a much smaller data structure, but what the heck
    for (unsigned int k=0; k<oligoGraph.size(); k++){
        oligoGraph.at(k).ED.resize(oligoGraph.size());
         for (unsigned int x=0; x<oligoGraph.size(); x++){
            string temp = oligoGraph.at(k).composition;
            pDist tempDist;
            temp=oligoGraph.at(x).composition;
                tempDist.dist =editDist(oligoGraph.at(k).composition, oligoGraph.at(x).composition);
                tempDist.numNode=x;
                if ( tempDist.dist >=MINEDIT)
                    oligoGraph.at(k).ED.push_back(tempDist);
        }
    }
*/

//Finding the cliques in a graph is NP-hard. I believe fixed  graphs is not, but hey!


//For nNow we force brute and test all our groups.


//For every oligo

 /* for (unsigned int k=0; k< oligoGraph.size();k++){

      //Generate every possible combination.
        cout <<oligoGraph.at(k).composition << "\t";

    } */



//pGraphs tempGraph;
/*
for (int i =0; i < tempCliques.size(); i++){
    tempGraph= tempCliques.at(i);
     cout << "Clique " << i << " ";
    for( int k =0 ;k< tempGraph.nodes.size(); k++ ){



        cout <<  oligoGraph.at(tempGraph.nodes.at(k)).composition  << " " ;

    }
    cout << endl;

}
*/


/*
  for (unsigned int k=0; k< oligoGraph.size();k++){
        cout <<oligoGraph.at(k).composition << "\t";
    }



cout <<endl;

  for (unsigned int k=0; k<oligoGraph.size(); k++){
        cout << oligoGraph.at(k).composition;
         for (unsigned int x=0; x<oligoGraph.size(); x++){

            cout << "\t" << oligoGraph.at(k).ED.at(x);
        }
        cout << endl;
    }

*/


    return 0;
}
//Assume same lenght string
int editDist( string curString, string otherString){

int dist=0;

    for (unsigned int i=0; i< curString.size(); i++){

        if(curString.at(i)!=otherString.at(i))
            dist++;
    }

return dist;
}



//Generate every permutation of lenght Index
//By reference if we need to optimize

vector<string> generateNucleoSubgroups(int index, vector<string> nucleo){

vector<string> temp_nucleo;

//End case
if (index ==1)
{
   temp_nucleo.push_back("A");
   temp_nucleo.push_back("G");
   temp_nucleo.push_back("C");
   temp_nucleo.push_back("T");
    return temp_nucleo;
}

//Recursion
    index--;
    nucleo = generateNucleoSubgroups(index, nucleo);

//generate substrings
    for( unsigned int i=0; i<nucleo.size(); i++){
          //Add A G C T
        temp_nucleo.push_back("A"+ nucleo.at(i));
        temp_nucleo.push_back("G"+ nucleo.at(i));
        temp_nucleo.push_back("C"+ nucleo.at(i));
        temp_nucleo.push_back("T"+ nucleo.at(i));
    }

return temp_nucleo;
}

pGraphs generateNucleoCliqueBetter(int index, const vector<oligoNode> oligoNodes, pGraphs oligoCliques , const int dist ){

    //Every Clique
    vector<pGraphs> temp_V_Clique;
    //One Clique
    pGraphs tempClique,bestClique;

    if (index==0){
        totalCount++;

        for (int k=0; k< oligoCliques.nodes.size();k++){
            cout << oligoNodes.at(oligoCliques.nodes.at(k)).composition << "\t";
    }
        cout <<endl;
        int tempCount = oligoCliques.nodes.size();
        //cout  << <<endl;
        //cout << "Found!" <<totalCount <<endl;
     if( (totalCount%5000) ==0)
        cerr << "Done " << totalCount <<"index"<<endl;

         return oligoCliques;
    }

int curNode, curCount, tempCount;

  if (oligoCliques.nodes.size()==0)
        curNode=0;
    else
        curNode= (oligoCliques.nodes.back()+1);
    index--;


for (int i =curNode; i<oligoNodes.size(); i++ ){
   tempClique= oligoCliques;
    //If this new add would not be a clique, we do not test it.
   if (isClique(&oligoCliques,oligoNodes, i, dist)){
            tempClique.nodes.push_back(i);
            tempClique=generateNucleoCliqueBetter(index,oligoNodes,tempClique, dist);
            tempCount= editCount(tempClique, oligoNodes);
            if ((tempCount> curCount) ||(bestClique.nodes.size()==0)){
                bestClique=tempClique;
                curCount= tempCount;
            }
   }
}

    return bestClique;

}



//We created every clique of one
//Then merge to every clique of two
//Etc
//Indexis the size of the cliques we want
//Dis is the min legal distance.
/*
vector<pGraphs> generateNucleoClique(int index, const vector<oligoNode> &oligoNodes, vector<pGraphs> oligoCliques , const int dist ){

    //Every Clique
    vector<pGraphs> temp_V_Clique;
    //One Clique
    pGraphs tempClique;


//Start case, make every clique of one
if (index ==1)
{
    for (int i=0; i<= oligoNodes.size(); i++){
       tempClique.nodes.clear();
       tempClique.nodes.push_back(i);

        temp_V_Clique.push_back(tempClique);
    }
    //return every clique of 1
    return temp_V_Clique;
}

//Recursion
    index--;
    //Get every legal clique of size index
    oligoCliques = generateNucleoClique(index, oligoNodes,  oligoCliques,dist);

//Add to every clique and validate.
    pGraphs *currentClique;
    for( unsigned int i=0; i<oligoCliques.size(); i++){

        //For every clique, we create a number of possible cliques and validate them.
        //We start at the last node, add every other possible node seperately and validate them.
        currentClique=&oligoCliques.at(i);

        for (int k=(currentClique->nodes.back()+1); k<oligoNodes.size(); k++ ){
            tempClique=*currentClique;
            //If this is a new legal clique, add it, otherwise do nothing
            if (isClique(currentClique,oligoNodes, k, dist ) )
            {
                tempClique.nodes.push_back(k);
                temp_V_Clique.push_back(tempClique);
                tempClique.nodes.clear();
            }
        }
    }

    cout << "Finished index" << index <<endl;
    cout << "There are " << temp_V_Clique.size() <<endl;


return temp_V_Clique;

}
*/

bool isClique(pGraphs *ourClique, const vector<oligoNode> &oligoNodes, int newNode , const int dist){

    int node1;
    int nodeDist;
    //Check distance between our newNode and every other clique
    for (int i=0; i<ourClique->nodes.size();i++ ){
        node1= ourClique->nodes.at(i);
        nodeDist= editDist(oligoNodes.at(node1).composition, oligoNodes.at(newNode).composition);
    if (nodeDist < dist)
        return false;
        //Does this new node connect to every other.
    }
return true;
}


//NUmber of edges with more then ED 2 in our Clique
int editCount(const pGraphs &ourClique,const vector<oligoNode> &oligoNodes){
  int node1, node2;
int count=0;

int temp= ourClique.nodes.size();

    for (int i=0; i< ourClique.nodes.size(); i++){
         node1= ourClique.nodes.at(i);
        for (int k=(i+1); k< ourClique.nodes.size(); k++){
            node2= ourClique.nodes.at(k);
            if ( editDist(oligoNodes.at(node1).composition, oligoNodes.at(node2).composition) > 2)
                count++;
        }
    }

    return count;
}

