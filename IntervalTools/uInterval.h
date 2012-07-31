#ifndef UINTERVAL_INCLUDED
#define UINTERVAL_INCLUDED

#include "uRegion.h"
#include "utility.h"
#include <string>

class uInterval : public uGenericNGS
{
    public:
        uInterval();
        uInterval(std::string chr, int start, int end,std::string extra);
        uInterval(uGenericNGS);

        void setExtra(std::string pExtra){m_Extra=pExtra;};
        std::string getExtra()const {return m_Extra;}  ;
        void writeAsTabFile(std::ostream& stream);

    protected:

    // User friendly name of our region.

    private:
         std::string m_Extra;
};

class uIntervalChrom : public uGenericNGSChrom<uInterval>
{
public :

    void writeAsTabFile(std::ostream& stream);


};

class uIntervalExperiment: public uGenericNGSExperiment<uIntervalChrom, uInterval>{
     public:

    void loadFromTabFile(std::ifstream& stream);
    void writeAsTabFile(std::ostream& stream);
};


namespace factory
{
static uInterval makeIntervalfromTabString(const std::string tabString)
{
    {
        utility::Tokenizer tabLine(tabString);
        std::string chrm,extra;
        int start=0, end=0;
        try
        {
            tabLine.NextToken();
            chrm = tabLine.GetToken();
            tabLine.NextToken();
            start = utility::stringToInt(tabLine.GetToken());
            tabLine.NextToken();
            end = utility::stringToInt(tabLine.GetToken());
            while(tabLine.NextToken())
                extra+= tabLine.GetToken();
            return uInterval(chrm,start,end,extra);
        }
        catch(...)
        {
            throw;
        }
    }
}

}

#endif // UINTERVAL_INCLUDED
