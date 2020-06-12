
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>

#include "random.h"
#include "functions.h"
using namespace std;


int main (int argc, char *argv[]){
//Generatore di numeri casuali
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   rnd.SaveSeed();
   
   
   /*****************************************************/
  double  L=1.0;
  double  d=4.0;
  int     M=100000; //numero lanci tot
  int     N=100; //numero di blocchi
  double  S=M/N;//numero di lanci per blocco
  int     N_h=0; //numero lanci positivi
  
  vector<int> x(N);
  iota(x.begin(), x.end(), 1);
  vector<double> AV(N, 0.0);
  vector<double> AV2(N, 0.0);
  vector<double> sum_prog(N, 0.0);
  vector<double> su2_prog(N, 0.0);
  vector<double> err_prog(N, 0.0);
  
  vector<double> x0; // vettore dei centri
  vector<double> xt;
  vector<double> yt; //xt e yt vettori per calcolare theta
  vector<double> t; //vettore theta

  
  for(int i=0; i<M; i++){
    x0.push_back(rnd.Rannyu(0, d));
  }
  
  for (int i=0; i<M; i++){ //estrazioni punti uniformi nella circonferenza unitaria
    double x1;
    double y1;
    do{
      x1=rnd.Rannyu(-3,3); //la scelta di -3,3 è stata fatta per avere uniformità dentro la circonferenza
      y1=rnd.Rannyu(-3,3);
    }while(pow(x1,2)+pow(y1,2)>1);
    xt.push_back(x1);
    yt.push_back(y1);
  }
  for (int i=0; i<M; i++){
    t.push_back(atan(yt[i]/xt[i])); //carico il vettore dei theta
  }
  
  
  for (int i=0; i<N; i++){
    int N_h=0; //lanci successo
    int k=i*S;
    N_h=Count(d, L, S, k, x0, t);
    AV[i]=double(N_h)/S;
    AV2[i]=pow(AV[i], 2);
  }
  
  for (int i=0; i<N; i++){
    for(int j=0; j<i+1; j++){
      sum_prog[i]+= AV[j];
      su2_prog[i]+= AV2[j];
    }
    sum_prog[i]/=(i+1);
    su2_prog[i]/=(i+1);
    err_prog[i]=Error(sum_prog, su2_prog, i );
  }
  
  ofstream Output;
  Output.open("pi.out");
  if (Output.is_open()){
    for (int i=0; i<N; i++){
      Output << x[i]*S << "\t" << (2*L)/(sum_prog[i]*d)  << "\t" << (2*L)/(pow(sum_prog[i],2)*d)*err_prog[i] << endl;
      cout << (2*L)/(sum_prog[i]*d) << endl;
    }
  } else cerr << "PROBLEM: Unable to open pi.out" << endl;
  Output.close();
  
  return 0;
}

