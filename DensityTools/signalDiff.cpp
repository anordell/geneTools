#include <tclap/CmdLine.h>
#include "signalDiff.h"
#include "utility/utility.h"
#include <string>
#include <string.h>
#include <fstream>
#include <ostream>
#include "uTags.h"
#include "uRegion.h"

using namespace NGS;
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


    uParser samParser(&samTreatStream,"SAM");

    tagTreat.loadWithParser(samParser);
    tagTreat.sortSites();


    uParser samControlParser(&samControlStream,"SAM");
    tagControl.loadWithParser(samControlParser);
    tagControl.sortSites();

    uRegionExperiment regionExpCtrl,regionExpTreat;

    {uParser regParser(&regionAStream,"BED");

    regionExpCtrl.loadWithParser(regParser);
    }
    uParser regParser(&regionAStream,"BED");
    regionExpTreat.loadWithParser(regParser);

    regionExpCtrl.measureDensityOverlap(tagControl);

    regionExpCtrl.generateSignal(tagControl);

    regionExpTreat.measureDensityOverlap(tagTreat);
    regionExpTreat.generateSignal(tagTreat);

    generateDistanceScores(regionExpTreat,regionExpCtrl);
    //  ofstream Ctrl("Ctrl.txt");
    ofstream Treat("Treat.txt");
    // regionExpCtrl.writeAll(Ctrl);

    uWriter regWriter(&Treat,"BED3");
    regionExpTreat.writeWithWriter(regWriter);
}

/** \brief Measure various distances between same elements of the exp and adds them to the elements
 *
 * \param regionExpA uRegionExperiment& Data A
 * \param regionExpB uRegionExperiment& Data B with sites positions ordered identical too A
 * \return void
 *
 */
void generateDistanceScores( uRegionExperiment & regionExpA, uRegionExperiment & regionExpB)
{

    /**< Careful, this is temporary */
    float maxAll=0.0;

    for (auto chromitA= regionExpA.begin(); chromitA!=regionExpA.end() ; chromitA++)
    {

      /*  for (int i=0; i<chromitA->second.count(); i++)
        {

        } */



       for (auto seqIT = chromitA->second.begin();seqIT!=chromitA->second.end();seqIT++ ){
            auto signal =seqIT->getSignal();
            maxAll=std::max(maxAll,*max_element(signal.begin(),signal.end()));
       }

    }
    /**< Compare each pair of regions and generate it's distances */
    for (auto chromitA= regionExpA.begin(); chromitA!=regionExpA.end() ; chromitA++)
    {
        auto chromitB= regionExpB.getpChrom(chromitA->first);

         auto seqITB =   chromitB->begin();
         for (auto seqITA = chromitA->second.begin();seqITA!=chromitA->second.end();seqITA++ ){

            vector<float> signalA= seqITA->getSignal();
            vector<float> signalB = seqITB->getSignal();

            int diffCount= (seqITA->getCount()-seqITB->getCount());
            seqITA->setScore(diffCount,0);

            float dist= clustering::align_distance(seqITA->getSignal(), seqITB->getSignal(),maxAll,(seqITA->getLength()/10));
            seqITA->setScore(dist,1);
            float hausdorffScore = clustering::hausdorffTwoRegions(signalA,signalB);
            seqITA->setScore(hausdorffScore,2);

            float normeuclidean = clustering::norm_euclidean_dist(signalA,signalB);
            seqITA->setScore(normeuclidean,3);

            seqITB++;
         }
    /*    for (int i=0; i<chromitA->second.count(); i++)
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

            float dist= clustering::align_distance(elemA->getSignal(), elemB->getSignal(),maxAll,(elemA->getLength()/10));
            elemA->setScore(dist,1);
            float hausdorffScore = clustering::hausdorffTwoRegions(signalA,signalB);
            elemA->setScore(hausdorffScore,2);

            float normeuclidean = clustering::norm_euclidean_dist(signalA,signalB);
            elemA->setScore(normeuclidean,3);
        }*/

    }

}



