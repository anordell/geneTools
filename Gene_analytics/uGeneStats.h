#ifndef UGENESTATS_H_INCLUDED
#define UGENESTATS_H_INCLUDED

#include "uChipSeq.h"
#include "uGeneRef.h"
#include "uGene.h"
class uGeneStats{


    private:

    uGeneRef ourGeneRef;
    uChipSeq ourChipSeq;

    public:

    uGeneStats(uGeneRef ourRef, uChipSeq ourChip);
      uGeneStats();
     ~uGeneStats();

    void loadData(uGeneRef ourRef, uChipSeq ourChip);

    void printStats(std::ostream& out);
    //Count how central peaks are within size distance of the tss
    int countTss(int size);
    //Count how many a proximal promoter of specified size
    int countinProximalProm(int size);
    //Count how many central peak-points are in a gene
    int countInGene();
    //Count how many central peak-points are in an Exon
    int countInExon();
    //Count how many peaks overlap the TSS
    int countOverlapTss();
    //Count how many peaks overlap a Gene
    int countOverlapgene();
    //Return a list of genes near each peak within distance. 0 means unlimited.
    vector<uGene> getNearestGenes(int distance);

    vector<int> uGeneStats::getTssProximityWindow();

};




#endif // UGENESTATS_H_INCLUDED
