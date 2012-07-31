#ifndef GENERATEEPIALPHA_H_INCLUDED
#define GENERATEEPIALPHA_H_INCLUDED

#include <map>
#include <vector>
#include <string>

class uAlphabetChrom;
class uAlphabetExperiment;

void writeMap(std::map<std::string,std::vector<int>> genomeMap, std::ostream & out);
void generateEpiAlpha(int argc, char* argv[]);
void exploreAlphabet(std::map<int,char>& alphabet,  std::vector<int>::const_iterator itrCurBasic,  std::vector<int>::const_iterator itrCurEnd, char & curChar, int curValue);
std::vector<int> setID( std::vector<uAlphabetExperiment> & histList);
std::map<int,char> createAlphabet(const std::vector<int> & idList);
std::vector<uAlphabetExperiment> loadData(std::vector<std::string> , std::vector<std::string> ,std::map<std::string,int> );
std::map<std::string, std::vector<int>> overlapGenomeWithMarks( std::vector<uAlphabetExperiment> p_markData,std::map<int,char>,std::string);

#endif // GENERATEEPIALPHA_H_INCLUDED
