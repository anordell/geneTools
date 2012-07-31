#include <tclap/CmdLine.h>
#include "signalDiff.h"
#include "utility.h"
#include <string>
#include <string.h>
#include <fstream>
#include <ostream>
#include "uTags.h"
#include "uRegion.h"


using namespace std;

void signalDiff(int argc, char* argv[])
{
    TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
    /**< Declare and add arguments */
    TCLAP::UnlabeledValueArg<string>  command("command", "Possible values are : peakStat", false," ","String");
    TCLAP::ValueArg<std::string> regionAarg("i","interval","File containing genomic interval regions",true,"null","string",cmd);
    TCLAP::ValueArg<std::string> samArgTreat("t","treat","Treatment single end sam files path",true,"null","string",cmd);
    TCLAP::ValueArg<std::string> samArgControl("c","control","No Treatment single end sam files path",true,"null","string",cmd);


    cmd.add(command);
    /**< Parse */
    cmd.parse( argc, argv );
    /**< Assign */
    std::string regionAPath = regionAarg.getValue();
    std::string samTreatPath = samArgTreat.getValue();
    std::string samControlPath = samArgControl.getValue();

    ifstream samTreatStream;
    utility::loadStream(samTreatPath, samTreatStream);
    ifstream samControlStream;
    utility::loadStream(samControlPath, samControlStream);
    ifstream regionAStream;
    utility::loadStream(regionAPath, regionAStream);
    cerr << " Loading" <<endl;
    /**< Overlap and get density */
    uTagsExperiment tagTreat, tagControl;
    tagTreat.loadFromSam(samTreatStream);
    tagTreat.sortData();

    tagControl.loadFromSam(samControlStream);
    tagControl.sortData();

    uRegionExperiment regionExpCtrl,regionExpTreat;
    regionExpCtrl.loadFromTabFile(regionAStream);
    regionExpTreat=regionExpCtrl;

    regionExpCtrl.measureDensityOverlap(tagControl);

    regionExpCtrl.generateSignal(tagControl);

    regionExpTreat.measureDensityOverlap(tagTreat);
    regionExpTreat.generateSignal(tagTreat);

    generateDistanceScores(regionExpTreat,regionExpCtrl);
    //  ofstream Ctrl("Ctrl.txt");
    ofstream Treat("Treat.txt");
    // regionExpCtrl.writeAll(Ctrl);
    regionExpTreat.writeAll(Treat);
}

/** \brief Measure various distances between same elements of the exp
 *
 * \param regionExpA uRegionExperiment&
 * \param regionExpB uRegionExperiment&
 * \return void
 *
 */
void generateDistanceScores( uRegionExperiment & regionExpA, uRegionExperiment & regionExpB)
{

    /**< Careful, this is temporary */
    float maxAll=0.0;

    for (auto chromitA= regionExpA.first(); chromitA!=regionExpA.last() ; chromitA++)
    {

        for (int i=0; i<chromitA->second.count(); i++)
        {
            auto elemA=chromitA->second.getPSite(i);
            auto signal =elemA->getSignal();
            maxAll=std::max(maxAll,*max_element(signal.begin(),signal.end()));
        }

    }

    for (auto chromitA= regionExpA.first(); chromitA!=regionExpA.last() ; chromitA++)
    {
        auto chromitB= regionExpB.getpChrom(chromitA->first);

        for (int i=0; i<chromitA->second.count(); i++)
        {
            auto elemA=chromitA->second.getPSite(i);
            auto elemB=chromitB->getPSite(i);

            vector<float> signalA= elemA->getSignal();
            vector<float> signalB = elemB->getSignal();
            //  float maxA, maxB;
            //  maxA=*std::max_element(signalA.begin(), signalA.end());
            // maxB=*std::max_element(signalB.begin(), signalB.end());

            int diffCount= (elemA->getCount()-elemB->getCount());
            elemA->setScore(diffCount,0);

            float dist= clustering::align_distance(elemA->getSignal(), elemB->getSignal(),maxAll,(elemA->getLenght()/10));
            elemA->setScore(dist,1);
            float hausdorffScore = clustering::hausdorffTwoRegions(signalA,signalB);
            elemA->setScore(hausdorffScore,2);

            float normeuclidean = clustering::norm_euclidean_dist(signalA,signalB);
            elemA->setScore(normeuclidean,3);
        }

    }

}



