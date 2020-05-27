#include "random.h"
#include "functions.h"

using namespace std;


int main (int argc, char *argv[]){

  PRNGen G; //setting generatore
  Random * rnd=G.SetGen("Primes", "seed.in");
  
  double t=0;
  double S_0=100;
  double T=1;
  double K=100;
  double r=0.1;
  double sgm=0.25;
  double M=1E5;
  double N=100;
  double L=M/N;
  
  vector<double> AVC(N, 0.0);
  vector<double> AVC2(N, 0.0);
  vector<double> AVP(N, 0.0);
  vector<double> AVP2(N, 0.0);
  
  //calcolo caso diretto
  
  for(int i=0; i<N; i++){
    double sumC=0;
    double sumP=0;
    for(int j=0; j<L; j++){
      double S_T=S_direct(S_0, r, sgm, T, rnd);
      sumC+=C00(S_T, r, T, K);
      sumP+=P00(S_T, r, T, K);
      }
    AVC[i]=sumC/L;
    AVC2[i]=pow(sumC/L, 2);
    AVP[i]=sumP/L;
    AVP2[i]=pow(sumP/L, 2);
  }
  Eval_Print("C_direct.out", AVC, AVC2, N, L);
  Eval_Print("P_direct.out", AVP, AVP2, N, L);
  
  //calcolo caso discreto
  
  for(int i=0; i<N; i++){
    double sumC=0;
    double sumP=0;
    for(int j=0; j<L; j++){
      double S_T=S_discret(S_0, r, sgm, T, rnd, 100);
      sumC+=C00(S_T, r, T, K);
      sumP+=P00(S_T, r, T, K);
    }
    AVC[i]=sumC/L;
    AVC2[i]=pow(sumC/L, 2);
    AVP[i]=sumP/L;
    AVP2[i]=pow(sumP/L, 2);
  }
  Eval_Print("C_discretized.out", AVC, AVC2, N, L);
  Eval_Print("P_discretized.out", AVP, AVP2, N, L);
  
  
  
  return 0;
  }

