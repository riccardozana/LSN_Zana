#include "random.h"
#include "functions.h"

using namespace std;


int main (int argc, char *argv[]){

  PRNGen G; //setting generatore
  Random * rnd=G.SetGen("Primes", "seed.in");
  
  int M=1E4; //definisco variabili dell'esperimento
  int N=100;
  int L=M/N;
  int dim=3;
  int steps=100;
  
  vector<vector<double> > matrix(M); // matrice per registrare coordinate
  for (int i=0; i<M; i++)
    matrix[i].resize(dim, 0.0);
  vector<double> AV(N, 0.0); //vettori per calcolo media a blocchi
  vector<double> AV2(N, 0.0);
  vector<vector<double> > output_d(steps); //matrice risultati
  
  for(int s=0; s<steps; s++){
  
    for(int i=0; i<M; i++){ //carico la matrice con i passi
      int dir=int(rnd->Rannyu(0, 3));
      int pos;
      if(rnd->Rannyu(-1, 1)>0) pos=+1;
      else pos=-1;
      matrix[i][dir]+=pos;
    }
    
    for (int i=0; i<N; i++){ // calcolo l'incertezza ad ogni passo con il metodo a blocchi
      double sum=0;
      for (int j=0; j<L; j++){
        int k=j+i*L;
        sum+=R2(matrix[k]);
      }
      AV[i]=sum/L;
      AV2[i]=pow(sum/L, 2);
    }
    Eval_Store(output_d[s], AV, AV2, N, L);
  }
  
  ofstream Output;
  Output.open("discreto.out");
  if (Output.is_open()){
    for(int i=0; i<output_d.size(); i++){
      Output << i+1 << "\t";
      for(int j=0; j<output_d[i].size(); j++){
        Output << output_d[i][j] << "\t";
      }
      Output << endl;
    }
  } else cerr << "PROBLEM: Unable to open -> " << "discreto.out" << endl;
  Output.close();
  
  /*****************************************************/
  //Caso Continuo
  
  for(int i=0; i<matrix.size(); i++){ //azzero la matrice
    for(int j=0; j<matrix[i].size(); j++){
      matrix[i][j]=0.0;
    }
  }
  
  vector<vector<double> > output_c(steps); //vettore risultati continuo
  
  for(int s=0; s<steps; s++){
    
    for(int i=0; i<M; i++){ //carico la matrice con i passi
      double theta=rnd->ThetaUnif();
      double phi=rnd->Rannyu(0, 2*M_PI);
      matrix[i][0]+=sin(theta)*cos(phi);
      matrix[i][1]+=sin(theta)*sin(phi);
      matrix[i][2]+=cos(theta);
    }
    
    for (int i=0; i<N; i++){ // calcolo l'incertezza ad ogni passo con il metodo a blocchi
      double sum=0;
      for (int j=0; j<L; j++){
        int k=j+i*L;
        sum+=R2(matrix[k]);
      }
      AV[i]=sum/L;
      AV2[i]=pow(sum/L, 2);
    }
    Eval_Store(output_c[s], AV, AV2, N, L);
  }
  

  Output.open("continuo.out");
  if (Output.is_open()){
    for(int i=0; i<output_c.size(); i++){
      Output << i+1 << "\t";
      for(int j=0; j<output_c[i].size(); j++){
        Output << output_c[i][j] << "\t";
      }
      Output << endl;
      }
      } else cerr << "PROBLEM: Unable to open -> " << "continuo.out" << endl;
      Output.close();
  
  return 0;
  }

