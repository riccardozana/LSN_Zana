
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
  //DADO UNIFORME
  int M=10000;
  
  vector<double> Sn1(M,0.0);
  vector<double> Sn2(M,0.0);
  vector<double> Sn10(M,0.0);
  vector<double> Sn100(M,0.0);
  
  
  //Carico i vettori con N numeri random uniformi per cella
  for(int i=0; i<1; i++){
    for (int j=0; j<M; j++){
      Sn1[j]+=(int (rnd.Rannyu(0.0, 6.0))+1);
    }
  }
  for(int i=0; i<2; i++){
    for (int j=0; j<M; j++){
      Sn2[j]+=(int (rnd.Rannyu(0.0, 6.0))+1);
    }
  }
  for(int i=0; i<10; i++){
    for (int j=0; j<M; j++){
      Sn10[j]+=(int (rnd.Rannyu(0.0, 6.0))+1);
    }
  }
  for(int i=0; i<100; i++){
    for (int j=0; j<M; j++){
      Sn100[j]+=(int (rnd.Rannyu(0.0, 6.0))+1);
    }
  }
  
  //Divido i contenuti dei vettori per N

  for(int j=0; j<M; j++){
    Sn1[j]/=1.;
    Sn2[j]/=2.;
    Sn10[j]/=10.;
    Sn100[j]/=100.;
  }
  
  //Infine stampo i valori
  
  ofstream Output;
  Output.open("unif_dice.out");
  if (Output.is_open()){
    for (int i=0; i<M; i++){
      Output << Sn1[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn2[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn10[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn100[i] << endl;
    }
  } else cerr << "PROBLEM: Unable to open unif_dice.out" << endl;
  Output.close();
  
  /*******************************************/
  //DADO ESPONENZIALE
  
  for(int j=0; j<M; j++){ //azzero i vettori
    Sn1[j]=0.0;
    Sn2[j]=0.0;
    Sn10[j]=0.0;
    Sn100[j]=0.0;
  }
  
  //Carico i vettori con N numeri random uniformi per cella
 
  for(int i=0; i<1; i++){
    for (int j=0; j<M; j++){
      Sn1[j]+=(rnd.Expo(1));
    }
  }
  for(int i=0; i<2; i++){
    for (int j=0; j<M; j++){
      Sn2[j]+=(rnd.Expo(1));
    }
  }
  for(int i=0; i<10; i++){
    for (int j=0; j<M; j++){
      Sn10[j]+=(rnd.Expo(1));
    }
  }
  for(int i=0; i<100; i++){
    for (int j=0; j<M; j++){
      Sn100[j]+=(rnd.Expo(1));
    }
  }
  
  //Divido i contenuti dei vettori per N
  
  for(int j=0; j<M; j++){
    Sn1[j]/=1.;
    Sn2[j]/=2.;
    Sn10[j]/=10.;
    Sn100[j]/=100.;
  }
  
  //Infine stampo i valori
  
  Output.open("expo_dice.out");
  if (Output.is_open()){
    for (int i=0; i<M; i++){
      Output << Sn1[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn2[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn10[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn100[i] << endl;
    }
  } else cerr << "PROBLEM: Unable to open expo_dice.out" << endl;
  Output.close();
  
  /*******************************************/
  //DADO LORENTZIANO
  
  for(int i=0; i<1; i++){
    for (int j=0; j<M; j++){
      Sn1[j]+=(rnd.Lorentz(0,1));
    }
  }
  for(int i=0; i<2; i++){
    for (int j=0; j<M; j++){
      Sn2[j]+=(rnd.Lorentz(0,1));
    }
  }
  for(int i=0; i<10; i++){
    for (int j=0; j<M; j++){
      Sn10[j]+=(rnd.Lorentz(0,1));
    }
  }
  for(int i=0; i<100; i++){
    for (int j=0; j<M; j++){
      Sn100[j]+=(rnd.Lorentz(0,1));
    }
  }
  
  //Divido i contenuti dei vettori per N
  
  for(int j=0; j<M; j++){
    Sn1[j]/=1.;
    Sn2[j]/=2.;
    Sn10[j]/=10.;
    Sn100[j]/=100.;
  }
  
  //Infine stampo i valori
  
  Output.open("lore_dice.out");
  if (Output.is_open()){
    for (int i=0; i<M; i++){
      Output << Sn1[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn2[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn10[i] << endl;
    }
    for (int i=0; i<M; i++){
      Output << Sn100[i] << endl;
    }
  } else cerr << "PROBLEM: Unable to open lore_dice.out" << endl;
  Output.close();
  
  
  
  
  
  
  
  
  
/*
  for (int i=0; i<M; i++){
    cout << Sn100[i] << endl;
  }
  */
  
  return 0;
}

