#include "scorenorm.h"
#include "functions.h"
//Establish a basic density mean and variance
//Then iterative score, normalize and decompose the files until every region left is below a given density/score threshold


using namespace std;

void scoreAndNormalize(int argc, char* argv[])
{


    if(argc <=3)
    {
        cerr<<"Program signature is -f <filepath> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsComplete=0]";
        return;
    }
    int threshold=0;
    bool isComplete=false;
    string pathname, chromName;
    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-f")==0)
        {
            // We know the next argument *should* be the filename:
            pathname = argv[i + 1];
        }
        else if (strcmp(argv[i],"-c")==0)
        {
            chromName = argv[i + 1];
        }
        else if (strcmp(argv[i],"-n")==0)
        {
            threshold = atoi(argv[i + 1]);
        }
        else if(strcmp(argv[i],"-p")==0)
        {
            int temp = atoi(argv[i+1]);
            if (temp)
                isComplete=true;
        }
    }
    //If we did not input file or chrom name
    if ((pathname.size()==0)||(chromName.size()==0))
    {
        cerr<<"Program signature is -f <filepath> -c <chromosomename> -n [Neighborhood threshold=0] -p [IsPe=0]";
        return;
    }

    //Finished Parsing parameters.
    uNucleoBin ourBin;
    vector<uRegion> ourRegions;
    //Load our Sam File
    {
        uTagsExperiment ourExp;
        ifstream inputStream;
        inputStream.open(pathname.c_str());
        //Load our Sam file
        if (inputStream.good())
            ourExp.loadSamHeader(inputStream);
        else
        {
            cerr << "Error loading in ScoreAndNormalize";
            abort();
        }
        //Loaded, generate our Density map per bp
        ourBin= loadDensityFromTags(ourExp, threshold,chromName, isComplete);
        cerr << "Finished density bins, Starting map Regions" <<endl;
        //Make our first map based on Nature score
        ourBin.generateSMap();
        cerr << "Finished map, starting Regions" <<endl;
        //Make Regions ( check this)
        ourRegions = makeRegionsfromBin(ourBin,ourExp);
        //Get SD and Mean
        statsStruct ourSD=returnSdfromRegions(ourRegions);
        cerr <<"Measuring Sum and Mean" << endl;
        cerr <<"SD is" <<ourSD.sd<< endl;
        cerr <<"Mean is" <<ourSD.mean<< endl;

        utility::pause_input();

        decompose(pathname, chromName, threshold, isComplete, ourSD);
    }
}

void decompose(string pathname,string chromName, int threshold, bool isComplete, statsStruct pSD)
{
        //From here, we run our decomp passes
        bool decomp=true;
        string currentpath=pathname;
        int passCount=0;
        while(decomp){

            //Load our tags, we will keep removing and using this.
          uTagsExperiment passData;
          cerr << "Loading following sam " << currentpath <<endl;
          ifstream currentData(currentpath);
          passData.loadFromSam(currentData);
          vector<int> depthVec;
          cerr <<"Starting pass " <<(passCount+1)<<endl;
            {
                //Set our IO for this pass
                ofstream ofLeft, ofDec, ofReg;
                ofDec.open(utility::concatStringInt("decomp.",passCount+1)+".sam");
                ofLeft.open(utility::concatStringInt("leftover.",passCount+1)+".sam");

                //Data to load
                auto ourRegions=getAndNormRegion(&passData,chromName,threshold,isComplete,pSD );
                writeRegNuclAsBedGraph(&ourRegions, (utility::concatStringInt("regions.",passCount)+".bedgraph"));

                ourRegions=filterRegions(ourRegions, 0.1f);
                //Why?
                vector<bedScores> ourBedScore=getBedFromEnlargedRegions(ourRegions, 30);
                //Find our Maxima.
                //Ok, but we could straight out use region middles
                ourBedScore=getMaximaFromBedscores(ourBedScore);
                //Take only the highest maxima from 150 bp regions and make a region out of it
                ourRegions=filterNucleoMaxima(&passData, chromName,ourBedScore);
                cerr << "There are " << ourRegions.size() << " regions after filtering for competitive max" <<endl;

                cerr << "Assignement depth map" <<endl;
                //Decompose by removing the 30 bp each side of every maxima from our Sam file.
                auto pChrom=passData.getpChrom(chromName);
                cerr << "Decomposing our maxima" <<endl;
                uTagsExperiment decompData=decompPass(pChrom,ourRegions);
                //Here, passe our regions and do +1 on iteration depth;
                cerr << "Sorting our decomposed data" <<endl;
                decompData.sortData();

                if (decompData.count()>1)
                {
                    passData.writeSamToOutput(ofLeft);
                    decompData.writeSamToOutput(ofDec);
                     //Next passs
                    currentpath=(utility::concatStringInt("leftover.",passCount+1)+".sam");
                }
                else
                    decomp=false;

                 passCount++;
            }
        }
}


//Eventually, change this to not be chrom dependent
 vector<uRegion>  getAndNormRegion(uTagsExperiment* pTags,string pChromname,  int pThreshold, bool pIsComplete, statsStruct pSD ){

    uTagsExperiment decompReturn;
    uNucleoBin curBin;

    //Etablish the number of tag centroids at each position
    curBin.setChromName(pChromname);
    curBin.generateDensityBin(pTags, pTags->getpChrom(pChromname)->getChromSize(), pThreshold,  pIsComplete);

    //Run the formulae from the Nature paper, giving us a specificity score at each position
    curBin.generateSMap();
    //Note that normally, our EXp here is not used..
    uTagsExperiment ourExp;
    auto ourRegions = makeRegionsfromBin(curBin,ourExp);
    //We normalize every score based on our SD and mean
    Normalize((pSD.sd/2), pSD.mean, ourRegions);
    return ourRegions;
}

void writeRegNuclAsBedGraph(vector<uRegion> * vecRegion, string ofRegPath){

            vector<uRegion>::iterator regit,ic;
            ofstream ofReg(ofRegPath.c_str());
            for ( regit=vecRegion->begin(); regit < vecRegion->end(); regit++ )
            {
                ofReg << regit->getChr() << "\t" << (regit->getStart()+29) << "\t" <<  (regit->getEnd()-30)<< "\t" << regit->getScore() << endl;
            }
            ofReg.close();
}
