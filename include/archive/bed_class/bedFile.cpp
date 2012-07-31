//Include various utility functions as well as the bedFile class.
//Might split in two in the future.


#include "bedFile.h"


using namespace std;


//int HAPPY=0;

bedFile::bedFile() {  // default constructor

    isWindowMapped=false;
    isIslandMapped=false;
    isTagMapped=false;

   for (int i=0; i<CHROMCOUNT;i++){
        genomeSumArray[i].wSize=-1;
        genomeSumArray[i].ChromNum=-1;
        genomeSumArray[i].ChromSize=-1;
    }
}

bedFile::~bedFile() {  // default constructor

}


void bedFile::reset(){

    isWindowMapped=false;
    isIslandMapped=false;
    isTagMapped=false;

    tagfile.clear();

    for (int i=0; i<CHROMCOUNT;i++){
        genomeSumArray[i].tList.clear();
        genomeSumArray[i].iList.clear();
        genomeSumArray[i].wList.clear();
        genomeSumArray[i].wSize=-1;
        genomeSumArray[i].ChromNum=-1;
        genomeSumArray[i].ChromSize=-1;
    }
}

ChromSum* bedFile::getChrom(int chromNumber) {

if (chromNumber>CHROMCOUNT)
    return 0;

return &genomeSumArray[chromNumber];
}

const Vec_Tag* bedFile::getTagVec(){

     return &tagfile;
 };

//Utility function for our tags


void bedFile::swapTags(Vec_Tag* tag_data,int posA, int posB){

    TagRead tempTag;
    tempTag=tag_data->at(posA);
    tag_data->at(posA)=tag_data->at(posB);
    tag_data->at(posB)=tempTag;

    return;
}





string bedFile::getChromfromint(int chromNum){


    string outputChrom;

    if (chromNum>=CHROMX_POS){
               if (chromNum== CHROMX_POS)
                    outputChrom="chrX";
                if (chromNum== CHROMY_POS )
                    outputChrom="chrY";
                if (chromNum== CHROMM_POS )
                    outputChrom="chrM";

           }
           else{
               outputChrom= "chr"+convertInt(chromNum+1);

           }
    return outputChrom;
}


int bedFile::getChromfromString(string chromName){

    int chromNum;
    //Ignore the first three letters
   // cout << chromName;
    chromName = chromName.substr(3);
    chromNum  = atoi(chromName.c_str());
    if ((chromNum==0)||(chromName.length()>4)){
        if (chromName=="M")
            chromNum=CHROMM_POS;
            else
            if (chromName=="X")
                chromNum=CHROMX_POS;
                else
                if (chromName=="Y")
                    chromNum=CHROMY_POS;
                    else
                    return -1;

     }
     else{
    chromNum--; }

    return chromNum;
}

//Wrapper so we don't aslways ask for the string
//User inputs string in this one
bool bedFile::openReadTags(){

    string str_filename;

    cout<<"Please enter file name containing alligned tags\n";
    cout<<"File must contain three columns, chrom, chromStart and chromEnd\n";
    cin >> str_filename;


 return readTags(str_filename);
}

//Need to add reading validation on this one at some point.
bool bedFile::readTags(string str_filename){

  ifstream i_bedfile;
TagRead tempRead;
  string  tempChrom;

    i_bedfile.open(str_filename.c_str() );
    if (!i_bedfile.is_open()){
          return false;
    }
    i_bedfile >> tempChrom;
    if (i_bedfile.eof())
            return false;
        else
        {
            i_bedfile.close();
            i_bedfile.open(str_filename.c_str() );
        }




  while (!i_bedfile.eof())
    {
        i_bedfile >> tempChrom;
         if( i_bedfile.eof() ) break;

        tempRead.chrom=getChromfromString(tempChrom);


            i_bedfile >> tempRead.chromStart;
            i_bedfile >> tempRead.chromEnd;

            //Ignore rest of line
            i_bedfile.ignore(5000, '\n');
         if (tempRead.chrom>=0){
             tagfile.push_back(tempRead);
        }

    }
    i_bedfile.close();


return true;
}

//Initialise our vectors to the window sizes

bool bedFile::initiWindowSpace(int wSize){
//For every Chromosome
      for (int i=0;i < CHROMCOUNT; i++)
    {
        //Easier to read the syntax if we assign the current chrom element to a point.
        ChromSum *pCurChrom;
        ChromWindow tempWindow;
        pCurChrom = &genomeSumArray[i];

        //Clear our data
        pCurChrom->wList.clear();

        //Basic window size
        pCurChrom->wSize=wSize;
        //Based on hg18 constants.
        pCurChrom->ChromSize=CHROMSIZE[i];
        pCurChrom->ChromNum=i;

        //Assign every window it's size, basic count and start position
        for(int k=0; k<( (pCurChrom->ChromSize)/wSize); k++){
            tempWindow.wSize= wSize;
            tempWindow.wCount=0;
            tempWindow.wStart=(k*wSize);
            pCurChrom->wList.push_back(tempWindow);


        }
        //If the chromosome does not neatly fit in our window size
        if ((pCurChrom->ChromSize%wSize)!=0)
        {
            tempWindow.wCount=0;
            //Final window equals the rest of our division.
            tempWindow.wSize=pCurChrom->ChromSize%wSize;
            tempWindow.wStart=pCurChrom->ChromSize/wSize;
            pCurChrom->wList.push_back(tempWindow);
        }
    }
    return true;
}

//All we want is to split our tags into there respective chromosomes, we do this.
bool bedFile:: mapTagstoChrom(){

//We map each tag into the chromosome list
 for( unsigned int i=0; i<tagfile.size();i++){
        int chromNumber;
        TagRead tempRead;
        ChromSum *pCurChrom;
        //What chromosome?
        tempRead=tagfile.at(i);
        chromNumber = tempRead.chrom;
        //Map pointer to that chromosome file.
        pCurChrom = &genomeSumArray[chromNumber];
        //add the tag to the list of tags on that chromosome.
        pCurChrom->tList.push_back(tempRead);
    }

    isTagMapped=true;
    return true;
}


//Map our tagfile into specific windows as determined by input.
//Cleans previous data.
bool bedFile::mapTagstoWindows(int wSize){

//If invalid parameters or data.
if ((wSize==0)||(tagfile.size()==0))
        return false;

//First we initialise our vector to window sizes.
    initiWindowSpace(wSize);

//Then we map each tag into one or several windows.
 for( unsigned int i=0; i<tagfile.size();i++){
        unsigned int chromNumber, windowPos;
        TagRead tempRead;
        ChromSum *pCurChrom;
        //What chromosome?
        tempRead=tagfile.at(i);
        chromNumber = tempRead.chrom;
        //Map pointer to that chromosome file.
        pCurChrom = &genomeSumArray[chromNumber];
        //add the tag to the list of tags on that chromosome.
        pCurChrom->tList.push_back(tempRead);

        windowPos = tempRead.chromStart / pCurChrom->wSize;

        //If the window start would be after our number of windows, problem!!!!
        if (windowPos >pCurChrom->wList.size() )
            return false;

    //Assign first tag
        pCurChrom->wList.at(windowPos).wCount++;


        //For now we map each tag to a single window.. for no particular reason since Cuddapah et al do it.
    //Does our Tag extend beyond our window?
    //We assume our reads map correctly.
        /*while(tempRead.chromEnd > (pCurChrom->wList.at(windowPos).wStart+pCurChrom->wList.at(windowPos).wSize)) {
            windowPos++;
            if (windowPos>=pCurChrom->wList.size())
                break;
            pCurChrom->wList.at(windowPos).wCount++;
            //In case of bug.
        }*/
    }
        isTagMapped=true;
        isWindowMapped = true;
    return true;
}



//We create islands bases on a tag threshold and a space of one "free" window.
bool bedFile::countIslands(const int minislanddensity, const int emptywindowLimit, const int minwindowdensity)
{
    unsigned int curTag;
    bool extend=true;
    int emptywindowcount=0;

    if (isWindowMapped==false)
        return false;




    for(int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        TagIsland tempIsland;
        pCurChrom = &genomeSumArray[i];
        //Clear old islands if there was any
        pCurChrom->iList.clear();

        curTag=0;
//God this is ugly. Refactor!!
        while (curTag < (pCurChrom->wList.size()) ) {
        //Start a potential  island
        extend=true;
        emptywindowcount=0;
            if (pCurChrom->wList.at(curTag).wCount>=minwindowdensity){
                tempIsland.iChrom=i;
                tempIsland.iCount=pCurChrom->wList.at(curTag).wCount;
                tempIsland.iStart=pCurChrom->wList.at(curTag).wStart;
                tempIsland.iStop=(tempIsland.iStart+pCurChrom->wList.at(curTag).wSize);
                tempIsland.iWindowCount=1;

                //Now we try to extend the window.
                while (extend==true)
                {
                    emptywindowcount=1;
                    curTag++;
                    if (curTag >=pCurChrom->wList.size() ){
                        extend=false;}
                    else
                    {
                        if (pCurChrom->wList.at(curTag).wCount>=minwindowdensity)
                        {
                            //If more then one tag in the window, we extend the island
                            tempIsland.iCount+=pCurChrom->wList.at(curTag).wCount;
                            tempIsland.iStop=(pCurChrom->wList.at(curTag).wStart+pCurChrom->wList.at(curTag).wSize);
                            tempIsland.iWindowCount++;
                        }
                        else
                        while(emptywindowcount<emptywindowLimit)
                           {
                                curTag++;
                                 if  (curTag >=pCurChrom->wList.size())
                                    {
                                    extend=false;
                                     emptywindowcount=emptywindowLimit;
                                    }
                                 //No we add and keep extending
                                 //( pCurChrom->wList.at(curTag).wCount<minwindowdensity )
                                else
                                {
                                    if (pCurChrom->wList.at(curTag).wCount>=minwindowdensity)
                                    {

                                        tempIsland.iCount+=pCurChrom->wList.at(curTag).wCount;
                                        tempIsland.iStop=(pCurChrom->wList.at(curTag).wStart+pCurChrom->wList.at(curTag).wSize);
                                       //Current window and the one we skipped.
                                        tempIsland.iWindowCount+=(emptywindowcount+1);
                                        emptywindowcount=emptywindowLimit;
                                    }
                                    else{
                                        emptywindowcount++;
                                        if (emptywindowcount==emptywindowLimit)
                                            extend=false;
                                        }
                                }
                           }
                    }
                }

            //Here we validate our window.
            //If more  then threshold tags validated
            if (tempIsland.iCount>minislanddensity )
                pCurChrom->iList.push_back(tempIsland);
            //Then go to next tag
            }

            curTag++;

        }//End when last tag checked

    //Change Chromosome
    }//End when last chromosome done

    isIslandMapped=true;
return true;
}

//Quisort implementation for sorting our chromosome.
//Still fairly long and could use some improvements.
void bedFile::sortChromTags(){

    ChromSum *pCurChrom;

      for(int i=0; i<CHROMCOUNT; i++){

        pCurChrom = &genomeSumArray[i];

        if (pCurChrom->tList.size()>1){
          QuickSortChrom(pCurChrom,0, (pCurChrom->tList.size()-1) );
          //QuickSortChrom(pCurChrom,12, (pCurChrom->tList.size()-1) );
        }
          //  QuickSortChrom(pCurChrom,0, (pCurChrom->tList.size()-1) );

     // cout <<"Finished Chrom "<<i<<"\n";
	}


}



 void bedFile::addRead(TagRead  tagtoadd){
     tagfile.push_back(tagtoadd);
return;
 }

void bedFile::addReadtoChrom(TagRead readtoAdd, int Chrom){

  try {


     if ( Chrom > CHROMCOUNT)
         throw (const char*)"Inserting into an illegal Chrom: Function AddReadtoChrom";


         genomeSumArray[Chrom].tList.push_back(readtoAdd);

  }
    catch (const char* msg) {
      cout << msg << endl;
      cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
      abort();
    }


return;
}

void bedFile::outputBedTags(string o_name){

    ofstream o_output;
    o_output.open((o_name.c_str()));
    string outputChrom;


    for (int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        pCurChrom = &genomeSumArray[i];

       for (unsigned int k=0; k <pCurChrom->tList.size();k++){

            outputChrom= getChromfromint(pCurChrom->tList.at(k).chrom);

            o_output << outputChrom << "\t";
            o_output << pCurChrom->tList.at(k).chromStart<< "\t";
            o_output << pCurChrom->tList.at(k).chromEnd << "\n";

            }

    }
}


void bedFile::outputWindows(string o_name){

    ofstream o_output;
    o_output.open((o_name.c_str()));
    string outputChrom;


    for (int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        pCurChrom = &genomeSumArray[i];

       for (unsigned int k=0; k <pCurChrom->wList.size();k++){

            outputChrom= getChromfromint(pCurChrom->ChromNum);

                o_output << outputChrom << "\t";
                o_output << pCurChrom->wList.at(k).wStart<< "\t";
                o_output << (pCurChrom->wList.at(k).wStart+pCurChrom->wList.at(k).wSize) << "\t";
                o_output << pCurChrom->wList.at(k).wCount<< "\n";
            }
    }
}


void bedFile::outputIslands(string o_name){

    ofstream o_output;
    o_output.open((o_name.c_str()));
    string outputChrom;

    for (int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        pCurChrom = &genomeSumArray[i];

       for (unsigned int k=0; k <pCurChrom->iList.size();k++){
           outputChrom= getChromfromint(pCurChrom->ChromNum);
                o_output << outputChrom << "\t";
                o_output << pCurChrom->iList.at(k).iStart<< "\t";
                o_output << pCurChrom->iList.at(k).iStop<< "\t";
                o_output << pCurChrom->iList.at(k).iWindowCount << "\t";
                o_output << pCurChrom->iList.at(k).iCount<<"\n";
            }
    }
}



void bedFile::QuickSortChrom(ChromSum* pCurChrom, int startIndex, int endIndex){

    TagRead pivot;

//HAPPY++;
//cout << HAPPY <<"\n";
   // cout << "Quicksort";

	               //pivot element first but should be random
	    int splitPoint;
	    if(endIndex > startIndex)                         //if they are equal, it means there is
	                                                      //only one element and quicksort's job
                                                    //here is finished
	    {
	        pivot = pCurChrom->tList.at(startIndex);

	        splitPoint = partition(pCurChrom, pivot.chromStart, startIndex, endIndex);
	                                                      //SplitArray() returns the position where
                                                        //pivot belongs to
	        pCurChrom->tList.at(splitPoint)=pivot;
	        QuickSortChrom(pCurChrom, startIndex, splitPoint-1);   //Quick sort first half
	        QuickSortChrom(pCurChrom, splitPoint+1, endIndex);    //Quick sort second half
	    }
}


//We steal someones partition function
int bedFile::partition(ChromSum* chromSum, int &pivot, int startIndex, int endIndex)
	{
	    int leftBoundary = startIndex;
	    int rightBoundary = endIndex;
        Vec_Tag* currentdata;



        currentdata= &chromSum->tList;
	    while(leftBoundary < rightBoundary)             //shuttle pivot until the boundaries meet
	    {
	         while( pivot < currentdata->at(rightBoundary).chromStart      //keep moving until a lesser element is found
	                && rightBoundary > leftBoundary)   //or until the leftBoundary is reached
	         {
	              rightBoundary--;                      //move left
	         }
	         swapTags(currentdata, leftBoundary, rightBoundary);
	         //PrintArray(array, ARRAY_SIZE);            //Uncomment this line for study

	         while( pivot >=currentdata->at(leftBoundary).chromStart        //keep moving until a greater or equal element is found
	                && leftBoundary < rightBoundary)   //or until the rightBoundary is reached
         {
              leftBoundary++;                        //move right
	         }
	         swapTags(currentdata,rightBoundary, leftBoundary);
	         //PrintArray(array, ARRAY_SIZE);            //Uncomment this line for study
	    }
	    return leftBoundary;                              //leftBoundary is the split point because
	                                                      //the above while loop exits only when
	                                                      //leftBoundary and rightBoundary are equal
	}
