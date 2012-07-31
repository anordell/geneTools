#include <string>
#include <vector>
#include <map>
#include "comparisonGroup.h"

/** @brief comparisonGroup
  *
  * @todo: document this function
  */

using namespace std;

/** @brief comparisonGroup
  *
  * @todo: document this function
  */
comparisonGroup::comparisonGroup(std::vector<std::string> geneIDVec, std::map<std::string, Presence> clusterMap)
{
    this->compMap=clusterMap;
    this->geneID=geneIDVec;
}

/** @brief comparisonGroup
  *
  * @todo: document this function
  */
comparisonGroup::comparisonGroup(std::vector<std::string> geneIDVec)
{
    this->geneID=geneIDVec;
}


comparisonGroup::comparisonGroup()
{
    //ctor
}

comparisonGroup::~comparisonGroup()
{

}



/** @brief getPresCount
  *
  * @todo: document this function
  */
int comparisonGroup::getPresCount(Presence presType)
{
    int counter=0;
    //map<std::string,int>::iterator it;

    for (auto itr = compMap.begin(); itr!=compMap.end(); itr++)
    {

        if (itr->second==presType)
            counter++;
    }


    return counter;
}

/** @brief returnCompResult
  *
  * @todo: document this function
  */
comparisonGroup comparisonGroup::returnCompResult(comparisonGroup pToCompare)
{

    vector<string> geneReturn=this->getGeneID();


    //Combine Gene names;
    vector<string> temp =pToCompare.getGeneID();
    geneReturn.insert( geneReturn.end(), temp.begin(), temp.end() );

    //Create with o9ur ID list
    comparisonGroup returnGroup(geneReturn);
//std::map<std::string, Presence>::iterator itr
    //For every Cluster in the new group, we compare and insert acceptable code.
    for(auto itr=pToCompare.begin(); itr!=pToCompare.end(); itr++)
    {
        {
            //If they do not match, set to MitchMach, else set to one of the two
            if (itr->second!=compMap[itr->first])
                returnGroup.addClusterData(itr->first,MIS);
            else
                returnGroup.addClusterData(itr->first,itr->second);
        }

    }

return returnGroup;
}

