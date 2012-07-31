#include <iostream>
#include <string>
#include <fstream>
#include  <vector>

using namespace std;


struct TagIsland{

 int iChrom;
 //NUmber of tags
 int iCount;
 int iStart;
 int iStop;
 int iWindowCount;
};


struct TagRead{

 //Only for alligned tags, not an actual bed format.
 int chrom;
 int chromStart;
 int chromEnd;
};

struct ChromWindow{

 //Chromosome
 int wStart;
 int wCount; //Number of aligned tags in the window.
 int wSize; //We keep for irregular steps.
};

struct ChromSum{

    vector<ChromWindow> wList; //Windows of the Chromosome
    vector<TagIsland>  wIslands;//Islands on the Chromose;
    int ChromNum; //Chromosome number -1 ( so Chr1= 0 ). See Constant def for CHROMM, X and Y
    int ChromSize;
    int wSize;

};


//Chrom Size
const int CHROM1 = 247249719;
const int CHROM2 = 242951149;
const int CHROM3 = 199501827;
const int CHROM4 = 191273063;
const int CHROM5 = 180857866;
const int CHROM6 = 170899992;
const int CHROM7 = 158821424;
const int CHROM8 = 146274826;
const int CHROM9 = 140273252;
const int CHROM10 = 135374737;
const int CHROM11 = 134452384;
const int CHROM12 = 132349534;
const int CHROM13 = 114142980;
const int CHROM14 = 106368585;
const int CHROM15 = 100338915;
const int CHROM16 = 88827254;
const int CHROM17 = 78774742;
const int CHROM18 = 76117153;
const int CHROM19 = 63811651;
const int CHROM20 = 62435964;
const int CHROM21 = 46944323;
const int CHROM22 = 49691432;
const int CHROMM  = 16571;
const int CHROMX  = 154913754;
const int CHROMY  = 57772954;


//Position in size vector
const int POSCHROMX = 22;
const int POSCHROMM = 23;
const int POSCHROMY = 24;

//There are 25 Chromosomes
const int CHROMCOUNT = 25;


const int TAGSPERWINDOW = 3;
//M, x then Y are are the last three


bool openReadTags(vector<TagRead>&,ifstream&,string );
bool readTags(vector<TagRead>& bedfile,ifstream& i_bedfile);
bool mapTags( vector<TagRead>,ChromSum[] );
bool countIslands( ChromSum genomeSumArray[]);

int main()
{



const int CHROMSIZE[]={CHROM1,CHROM2,CHROM3,CHROM4,CHROM5,CHROM6,CHROM7,CHROM8,CHROM9,CHROM10
                        ,CHROM11,CHROM12,CHROM13,CHROM14,CHROM15,CHROM16,CHROM17,CHROM18,CHROM19,CHROM20,CHROM21,
                        CHROM22,CHROMX,CHROMM,CHROMY};

    string str_filename, str_ofilename;
    vector<TagRead> bedfile;
  //  vector<int> chromSizes;
    ifstream i_bedfile;
    int windowSize;

    ChromSum genomeSumArray[25];


   if(!openReadTags(bedfile,i_bedfile, str_filename))
    {
        cout << "Program failed to open file. Press ENTER to exit";
        cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

        return 0;
    }

//Tag data is read, we create summary windows based on chrom sizes and choice of user

    cout << "Please input windows size. Suggested value between 100 and 2000\n";
    cin >> windowSize;


    for (int i=0;i < CHROMCOUNT; i++)
    {
        //Easier to read the syntax if we assign the current chrom element to a point.
        ChromSum *pCurChrom;
        ChromWindow tempWindow;

        pCurChrom = &genomeSumArray[i];
        pCurChrom->wSize=windowSize;
        pCurChrom->ChromSize= CHROMSIZE[i];
        pCurChrom->ChromNum=i;

        for(int k=0; k<( (pCurChrom->ChromSize)/windowSize); k++){
            tempWindow.wSize= windowSize;
            tempWindow.wCount=0;
            tempWindow.wStart=(k*windowSize);
            pCurChrom->wList.push_back(tempWindow);
        }
        //If the chromosome does not neatly fit in our window size
        if ((pCurChrom->ChromSize%windowSize)!=0)
        {
            tempWindow.wCount=0;
            //Final window equals the rest of our division.
            tempWindow.wSize=pCurChrom->ChromSize%windowSize;
            tempWindow.wStart=pCurChrom->ChromSize/windowSize;
            pCurChrom->wList.push_back(tempWindow);
        }

    }

    //We assign every tag to one or many windows. with 36 bp tags, this should only be two max.
    //We assumed our bedfiles are sorted correctly..otherwise insert BedSort here.

    if (mapTags(bedfile,genomeSumArray ) ==false){
      cout << "Program failed to read tags intro windows. Press ENTER to exit\n";
      cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );
        return 0;

    }


      cout << "SUCCESSS. Now we try to make islands\n";
      cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );



    countIslands(genomeSumArray);


    cout << "Islands made. Please input output file name. Warning, this will overwrite anything with the same name\n";
    cin >> str_ofilename;

    ofstream o_islands, o_stats;
    o_islands.open((str_ofilename.c_str()));
    str_ofilename=str_ofilename+".stats";
    o_stats.open((str_ofilename.c_str()));


    for (int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        pCurChrom = &genomeSumArray[i];

        o_stats << "Chromosome: "<<(i+1)<<" has "<<pCurChrom->wIslands.size()<<" islands.\n";
       for (int k=0; k <pCurChrom->wIslands.size();k++){


        //Output each Island;
        if (pCurChrom->wIslands.at(k).iWindowCount>=0){
        o_islands << pCurChrom->wIslands.at(k).iChrom << "\t";
        o_islands << pCurChrom->wIslands.at(k).iStart<< "\t";
        o_islands << pCurChrom->wIslands.at(k).iStop<< "\t";
        o_islands << pCurChrom->wIslands.at(k).iWindowCount << "\t";
        o_islands << pCurChrom->wIslands.at(k).iCount<<"\n";

        }
        }

    }



    cout << "Program stop. Press ENTER to exit";
    cin.ignore( std::numeric_limits<std::streamsize>::max(), '\n' );

return 0;
}


//Wrapper so we don't aslways ask for the string
//User inputs string in this one
bool openReadTags(vector<TagRead>& bedfile,ifstream& i_bedfile,  string str_filename){

    cout<<"Please enter file name containing alligned teags\n";
    cout<<"File must contain three columns, chrom, chromStart and chromEnd\n";
    cin >> str_filename;
    i_bedfile.open(str_filename.c_str() );
    if (!i_bedfile.is_open()){
          return false;
    }

 return readTags(bedfile, i_bedfile);
}

//Need to add reading validation on this one at some point.
bool readTags(vector<TagRead>& bedfile,ifstream& i_bedfile){

  int getChromfromString(string);

  TagRead tempRead;
  string  tempChrom;
  while (!i_bedfile.eof())
    {
        i_bedfile >> tempChrom;

        tempRead.chrom=getChromfromString(tempChrom);

        i_bedfile >> tempRead.chromStart;
        i_bedfile >> tempRead.chromEnd;

        //Ignore rest of line
        i_bedfile.ignore(255, '\n');
        bedfile.push_back(tempRead);

    }
return true;
}


int getChromfromString(string chromName){

    int chromNum;
    //Ignore the first three letters
    chromName = chromName.substr(3);
    chromNum  = atoi(chromName.c_str());
    if (chromNum==0){
        if (chromName=="M")
            chromNum=POSCHROMM;
        if (chromName=="X")
            chromNum=POSCHROMX;
        if (chromName=="Y")
            chromNum=POSCHROMY;
     }
     else{
    chromNum--; }

    return chromNum;
}


//Map tags into specified window sizes.
bool mapTags( vector<TagRead> bedfile, ChromSum genomeSumArray[] ){

 for( int i=0; i<bedfile.size();i++){
        int chromNumber, windowPos;
        TagRead tempRead;
        ChromSum *pCurChrom;
        //What chromosome?
        tempRead=bedfile.at(i);
        chromNumber = tempRead.chrom;
        //Map pointer to that file
        pCurChrom = &genomeSumArray[chromNumber];

        windowPos = tempRead.chromStart / pCurChrom->wSize;

        //If the window start would be after our number of windows, problem!!!!
        if (windowPos >pCurChrom->wList.size() )
            return false;

    //Assign first tag
        pCurChrom->wList.at(windowPos).wCount++;
    //Does our Tag extend beyond our window?
    //We assume our reads map correctly.
        while(tempRead.chromEnd > (pCurChrom->wList.at(windowPos).wStart+pCurChrom->wList.at(windowPos).wSize)) {
            windowPos++;
            if (windowPos>=pCurChrom->wList.size())
                break;
            pCurChrom->wList.at(windowPos).wCount++;
            //In case of bug.
        }
    }
    return true;
}

bool countIslands( ChromSum genomeSumArray[])
{
    int curTag;
    bool extend=true, emptywindowallowed=false;
    for(int i=0; i<CHROMCOUNT; i++){
        ChromSum *pCurChrom;
        TagIsland tempIsland;
        pCurChrom = &genomeSumArray[i];

        curTag=0;

        while (curTag < (pCurChrom->wList.size()) ) {
        //Start a potential  island
        extend=true;
        emptywindowallowed=false;
            if (pCurChrom->wList.at(curTag).wCount>=TAGSPERWINDOW){
                tempIsland.iChrom=i;
                tempIsland.iCount=pCurChrom->wList.at(curTag).wCount;
                tempIsland.iStart=pCurChrom->wList.at(curTag).wStart;
                tempIsland.iStop=(tempIsland.iStart+pCurChrom->wList.at(curTag).wSize);
                tempIsland.iWindowCount=1;

            //Now we try to extend the window.
                while (extend==true){
                    curTag++;
                    if (curTag >=pCurChrom->wList.size() ){
                        extend=false;}
                    else{
                        if (pCurChrom->wList.at(curTag).wCount>=TAGSPERWINDOW){
                            //If more then one tag in the window, we extend the island
                            tempIsland.iCount+=pCurChrom->wList.at(curTag).wCount;
                            tempIsland.iStop=(pCurChrom->wList.at(curTag).wStart+pCurChrom->wList.at(curTag).wSize);
                            tempIsland.iWindowCount++;
                        }
                        else{
                        //Is there only one empty window?
                            if (emptywindowallowed)
                            {
                                emptywindowallowed=false;
                                curTag++;
                                //Yes, we stop
                                 if ( (curTag >=pCurChrom->wList.size() )||( pCurChrom->wList.at(curTag).wCount<TAGSPERWINDOW )  )
                                 {extend=false;
                                 }
                                 //No we add and keep extending
                                else{
                                    tempIsland.iCount+=pCurChrom->wList.at(curTag).wCount;
                                    tempIsland.iStop=(pCurChrom->wList.at(curTag).wStart+pCurChrom->wList.at(curTag).wSize);
                                    tempIsland.iWindowCount++;
                                    }
                            }
                            else{


                                extend=false;
                            }
                        }
                    }
                }
                //Here we validate our window.
                //If fewer then 8 tags, invalide
                if (tempIsland.iCount>34 )
                    pCurChrom->wIslands.push_back(tempIsland);
            }
            //Then go to next tag
            curTag++;
        }


    }
return true;
}




