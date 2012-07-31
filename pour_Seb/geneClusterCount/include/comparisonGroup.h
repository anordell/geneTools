#ifndef COMPARISONGROUP_H
#define COMPARISONGROUP_H
#include <map>
enum Presence { ABS,PRES, MIS };
class comparisonGroup
{

  struct ClustData{
  std::string ClustID;
  Presence result;
  };
    public:
        comparisonGroup(std::vector<std::string> geneIDVec,std::map<std::string,Presence> clusterMap);
        comparisonGroup();
        comparisonGroup(std::vector<std::string> geneIDVec);
        virtual ~comparisonGroup();

      void addGeneID(std::string pgeneID){geneID.push_back(pgeneID);};
      void addClusterData(std::string ID, Presence result){ compMap.insert(std::pair<std::string, Presence>(ID, result));};
      void replaceClusterData(std::string ID, Presence result){compMap[ID]=result;};

      int getPresentCount(){return getPresCount(PRES);};
      int getAbsentCount(){return getPresCount(ABS);};
      int getMismatchCount(){return getPresCount(MIS);};
      int getGeneCount(){return geneID.size();};
      std::string getConcatGeneID(){
      std::string returnString;

        for (std::string x : geneID){

            returnString+=" ";
            returnString+=x;

            }
        return returnString;
      }
      std::vector<std::string> getGeneID(){return geneID;};


      std::map<std::string, Presence>::iterator begin(){return compMap.begin();};
      std::map<std::string, Presence>::iterator end(){return compMap.end();};
      comparisonGroup returnCompResult(const comparisonGroup pToCompare);

    protected:
    private:

        int getPresCount(Presence presType);
         /**< Contains the list of Gene ID we have included in this comparison */
        std::vector<std::string> geneID;
        /**< Map containing our comparisons, with string being the unique ID of the cluste and each entry the status of one cluster */
        std::map<std::string,Presence> compMap;

};
#endif // COMPARISONGROUP_H
