
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
   //Calcolo il valore medio
   int M=10000; //number of throws
   int N=100; //number of block
   int L=(M/N); //number of throws in each block
   
   vector<double> r;
      for (int i=0; i<M; i++){
         r.push_back(rnd.Rannyu());
      }//Uniform random number [0,1)
   
   vector<int> x(N);
   iota(x.begin(), x.end(), 1);
   vector<double> AV(N, 0.0);
   vector<double> AV2(N, 0.0);
   vector<double> sum_prog(N, 0.0);
   vector<double> su2_prog(N, 0.0);
   vector<double> err_prog(N, 0.0);
   
   for (int i=0; i<N; i++){
      double sum=0;
      for (int j=0; j<L; j++){
         int k=j+i*L;
         sum+=r[k];
      }
      AV[i]=sum/L;
      AV2[i]=pow(sum/L, 2);
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
   Output.open("mean.out");
   if (Output.is_open()){
      for (int i=0; i<N; i++){
      Output << x[i]*L << "\t" << sum_prog[i]  << "\t" << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open mean.out" << endl;
   Output.close();
   
   /*********************************+++++*/
   //Calcolo Dev Std
   
   for (int i=0; i<N; i++){
      sum_prog[i]=0;
      su2_prog[i]=0;
      err_prog[i]=0;
   }//azzero vettori
   
   for (int i=0; i<N; i++){
      double sum=0;
      for (int j=0; j<L; j++){
         int k=j+i*L;
         sum+=pow(r[k]-0.5,2);
      }
      AV[i]=sum/L;
      AV2[i]=pow(sum/L, 2);
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
   
   Output.open("dev.out");
   if (Output.is_open()){
      for (int i=0; i<N; i++){
         Output << x[i]*L << "\t" << sum_prog[i] -1.0/12.0  << "\t" << err_prog[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open dev.out" << endl;
   Output.close();
   
   /*********************************************/
   //Chi quadro Test
   
   int      n=10000;
            M=100;
   double   I=1./M;
   
   vector<int>counter (N,0);
   vector<double> sum_chi (N,0.0);
   
   for (int i=0; i<100; i++){
      r.clear(); //pulisco r
      int count=0;
      
      for (int j=0; j<n; j++){
         r.push_back(rnd.Rannyu()); //riempio r con numeri random
      }
      for(int j=0; j<n; j++){//conteggio le estrazioni per ogni bin
         while(r[j]>I*(count+1)){
            count++;
         }
         counter[count]++;//aumento di uno il conteggio nel bin "count"
         count=0; //riazzero il contatore
      }
      for(int j=0; j<M; j++){//registro il valore di chi quadro
         sum_chi[i]+=pow(counter[j]-n/M,2)/double(n/M);
         counter[j]=0;
      }
   }
   
   Output.open("chi2.out");
   if (Output.is_open()){
      for (int i=0; i<N; i++){
         Output << x[i] << "\t" << sum_chi[i] << endl;
      }
   } else cerr << "PROBLEM: Unable to open chi2.out" << endl;
   Output.close();
 

   return 0;
}

