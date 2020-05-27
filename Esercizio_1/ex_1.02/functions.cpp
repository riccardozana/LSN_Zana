#include "functions.h"


//Error function: statistical uncertainty estimation
double Error(std::vector<double> AV, std::vector<double> AV2, int n){
   if(n==0){
      return 0;
   }else{
      return sqrt((AV2[n]-pow(AV[n],2))/double(n));
   }
}


