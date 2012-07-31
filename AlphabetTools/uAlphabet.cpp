#include "uFormats.h"
#include "uAlphabet.h"


using namespace std;

/** \brief Default constructor
 */
uAlphabet::uAlphabet()
{
}

/** \brief Alphabet constructor, with start, end and chr
 *
 * \param ourchr std::string  : Chrom to set
 * \param ourstart int        : Start position
 * \param ourend int          : End position
 *
 */
uAlphabet::uAlphabet(std::string ourchr, int ourstart, int ourend)
{

    try {

        setChr(ourchr);
        setStartEnd(ourstart,ourend);
    }
    catch(...){
    cerr <<"Error in constructor of uAlphabet";
        throw;


    }

}

uAlphabet::uAlphabet(uGenericNGS tag)
{

    setChr(tag.getChr());
    setStartEnd(tag.getStart(),tag.getEnd());
}

void uAlphabetChrom::modifyAlphabet(vector<int> & ourAlphabet, int pID)
{

    applyOnAllSites([&](uAlphabet Tag)
    {

        int begin = Tag.getStart();
        int end = Tag.getEnd();

        for(int i=begin; i<=end; i++){


            if (i>=(int)ourAlphabet.size())
            {
                cerr <<i <<" from max of " <<ourAlphabet.size();
                utility::pause_input();
            }

            ourAlphabet.at(i)+=pID;


        }
    } );
}

void uAlphabetExperiment::modifyAlphabet(string curChrom, vector<int> & vecAlphabet)
{
        ExpMap[curChrom].modifyAlphabet(vecAlphabet, getID());

}


void uAlphabetExperiment::setID(int pID)
{
    m_IDa=pID;
    applyOnAllChroms([&](uAlphabetChrom mChrom)
    {

        mChrom.setID(pID);
    });


}


void uAlphabetExperiment::setExplain(std::string pExplain)
{

    m_Explanation=pExplain;
    applyOnAllChroms([&](uAlphabetChrom mChrom)
    {

        mChrom.setExplain(pExplain);

    });
}


void uAlphabetExperiment::loadFromTabFile(std::ifstream& stream)
{
    std::string tempString;
    while(!std::getline(stream, tempString).eof())
    {
        //  auto haha =  static_cast<_BASE_>(factory::makeNGSfromTabString(tempString));

        addSite(factory::makeAlphafromTabString(tempString));
    }

}


