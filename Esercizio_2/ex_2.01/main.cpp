

#include "random.h"
#include "functions.h"
#include "funzioneBase.h"
#include "integralMC.h"
using namespace std;


int main (int argc, char *argv[]){

  PRNGen G; //setto il generatore
  Random * rnd=G.SetGen("Primes", "seed.in");
  IntegralMC myMC(rnd);
  
  Cos * g= new Cos(M_PI/2, M_PI/2, 0.);

  int M=1E6;
  int N=100;
  int L=M/N;
  vector<double> AV(N, 0.0);
  vector<double> AV2(N, 0.0);
  
  for (int i=0; i<N; i++){
    AV[i]=myMC.IntegralAVE(0., 1., g, L);
    AV2[i]=pow(AV[i], 2);
  }
  
  Eval_Print("Iunif.out", AV, AV2, N, L); // Calcolo incertezza a blocchi e stampa
  
  /**************************************************/
  
  TaylorCos * t= new TaylorCos();
  Ratio * f= new Ratio(g, t);
  
  for (int i=0; i<N; i++){
    AV[i]=myMC.IntegralIMP(0., 1., t->Eval(0.), t, f, L);
    AV2[i]=pow(AV[i], 2);
  }
  
  Eval_Print("IImpSam.out", AV, AV2, N, L);
  
  return 0;
  }

