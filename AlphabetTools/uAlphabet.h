#ifndef UALPHABET_H
#define UALPHABET_H


#include "uRegion.h"

class uAlphabet : public uGenericNGS
{
    public:
        uAlphabet();
        uAlphabet(std::string chr, int start, int end);
        uAlphabet(uGenericNGS);
    protected:

    // User friendly name of our region.
};

class uAlphabetChrom : public uGenericNGSChrom<uAlphabet>
{

    public:

    void setID(int pID){m_IDa=pID;};
    int getID()const {return m_IDa;};

    void setExplain(std::string pExplain){m_Explanation=pExplain;};
    std::string getExplain()const {return m_Explanation;}  ;

    void modifyAlphabet(std::vector<int> & ourAlphabet, int pID);
    private:
    std::string m_Explanation;
    int m_IDa;
};

class uAlphabetExperiment: public uGenericNGSExperiment<uAlphabetChrom, uAlphabet>{
     public:

    void setID(int pID);
    void modifyAlphabet(std::string, std::vector<int> & );
    void setExplain(std::string pExplain);
    std::string getExplain(){return m_Explanation ;};

    int getID(){return m_IDa;};
    std::string m_Explanation;
    int m_IDa;

    void loadFromTabFile(std::ifstream& stream);

};




namespace factory
{
static uAlphabet makeAlphafromTabString(const std::string tabString)
{
    {

        utility::Tokenizer tabLine(tabString);

        std::string chrm;
        int start=0, end=0;
        try
        {
            tabLine.NextToken();
            chrm = tabLine.GetToken();
            tabLine.NextToken();
            start = utility::stringToInt(tabLine.GetToken());
            tabLine.NextToken();
            end = utility::stringToInt(tabLine.GetToken());

            return uAlphabet(chrm,start,end);
        }
        catch(...)
        {
            throw;
        }

    }


}

}


#endif
