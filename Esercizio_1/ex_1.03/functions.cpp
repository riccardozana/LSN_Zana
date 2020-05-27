#include "functions.h"


//Error function: statistical uncertainty estimation
double Error(std::vector<double> AV, std::vector<double> AV2, int n){
   if(n==0){
      return 0;
   }else{
      return sqrt((AV2[n]-pow(AV[n],2))/double(n));
   }
}

int Count (double d, double L, int M, int k, std :: vector<double> x0, std :: vector<double> t){
  int N_h=0;
  for(int i=k; i<(k+M); i++){
    if(x0[i]-abs(cos(t[i])*(L/2.))<0) N_h++;
    if(x0[i]+abs(cos(t[i])*(L/2.))>d) N_h++;
  }
  return N_h;
}

