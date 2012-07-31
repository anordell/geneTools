#ifndef CLUSTERING_H_INCLUDED
#define CLUSTERING_H_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
namespace clustering
{


//Distance between two set of points, using position in the vector as X
static inline float hauftsman(const std::vector<float> & vecA, const std::vector<float> & vecB)
{
    int sizeA = vecA.size(), sizeB = vecB.size();

    if (sizeA == 0 || sizeB == 0)
    {
        return 0.0;
    }

    float maxDistOverA= 0.0;

    int aPos=0, bPos=0;
    int norm= sizeA;
    for (auto valA= vecA.begin(); valA!=vecA.end(); valA++)
    {
//Todo : Formalize this
    bPos=0;
    float minDistOverB = 1000000;

  //  std::cout << "Just before "<< minDistOverB <<std::endl;
    for (const float& valB: vecB)
            {
                float diffx =(bPos-aPos);
               // std::cerr <<"X "<< aPos << " " << bPos <<std::endl;
               // std::cerr <<"Y "<< *valA << " " << valB <<std::endl;


                //if (diffx>0)
                 //   diffx=diffx/norm;
                diffx=pow(diffx,2);



                float diffy(valB-*valA );

                if (diffy!=0)
                    diffy=diffy*(norm/2);
                // std::cerr << "before abs "<< diffy << std::endl;
                diffy= pow(diffy,2);
                // std::cerr << "dist "<< diffx << " "<< diffy << std::endl;

                float dist = diffx + diffy ;
               // std::cerr << "dist "<< dist << std::endl;
                //if ((*valA!=valB)&&(bPos==aPos))
                //    pause_input();

                if (dist < minDistOverB)
                    minDistOverB = dist;
                bPos++;
            }
                // std::cout << "Assigning A to "<<minDistOverB <<std::endl;
            if (minDistOverB > maxDistOverA){

                maxDistOverA = minDistOverB;
            }
            aPos++;
    }

//Squared distance
    // std::cout << "Before return " <<std::endl;
    // std::cout <<maxDistOverA <<std::endl;
   // std::cout << "Previous val " <<std::endl;
    float returnval=0.0;
          //  maxDistOverA=maxDistOverA*100;
    if (maxDistOverA>0){
        returnval = sqrt(maxDistOverA);
       // returnval=returnval/100;
    }

    return returnval;
}

}

#endif // CLUSTERING_H_INCLUDED
