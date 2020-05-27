#include "random.h"
#include "functions.h"
#include "funzioneBase.h"
#include "integralMC.h"

using namespace std;


int main (int argc, char *argv[]){

  PRNGen G; //setting generatore
  Random * rnd=G.SetGen("Primes", "seed.in");
  
  vector<double> y(3, 0.0);
  vector<double> x(3, 0.0);
  
  
  
  /*****************************************************/
  //T-UNIF comportamento partenza lontano dall'origine
  /*****************************************************/

  
  //initial position
  y[0]=50;
  y[1]=50;
  y[2]=50;
  
  cout << "-------------------------" << endl;
  cout << "Initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  
  //grandezze simulazione
  double radius=1.55;
  double a0=1.; //unità riscalate con raggio bohr=1
  
  Psi_1_0_0 * psi1 = new Psi_1_0_0(a0);
  Retta * g=new Retta(1., 0.);
  
  //equilibrazione della posizione iniziale
  
  Equilibrate_unif_print (y, psi1, radius, 1000, rnd);
  cout << "Equilibrated initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  cout << "-------------------------" << endl;
  
  
  /*****************************************************/
  //T-UNIF
  /*****************************************************/
  
  //initial position
  y[0]=1/sqrt(3);
  y[1]=1/sqrt(3);
  y[2]=1/sqrt(3);
  
  cout << "-------------------------" << endl;
  cout << "Initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  
  //passi e blocchi simulazione
  int M=  1E6;
  int count=0;
  int N=  100;
  int L=  M/N;
  vector<double>   AV(N, 0.0);
  vector<double>  AV2(N, 0.0);
  
  //grandezze simulazione
  radius=1.55;
  a0=1.; //unità riscalate con raggio bohr=1
  
  //equilibrazione della posizione iniziale
  
  Equilibrate_unif (y, psi1, radius, 1000, rnd);
  cout << "Equilibrated initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  cout << "-------------------------" << endl;
  
  //calcolo dei valori della simulazione
  
  ofstream Output;
  Output.open("xyzpsi100.out");
  
  for(int i=0; i<N; i++){
    
    double sum=0;
    
    for(int j=0; j<L; j++){
      for(int k=0; k<y.size(); k++){
        x[k]=y[k];
      }
      
      T_unif(x, radius, rnd);
      double alpha=Eval_Alpha(y, x, psi1);
      Accept(alpha, y, x, count, rnd);
      sum+=g->Eval(sqrt(R2(y)));
      if(j%100==0)Output<< y[0] << "\t" << y[1] << "\t" << y[2] << endl;
      
    }
    
    AV[i] =sum/L;
    AV2[i]=pow(sum/L,2);
    
  }
  
  Output.close();
  Eval_Print("psi100.out", AV, AV2, N);
  
  
  Output.open("psi100_acceptance.out");
  Output << double(count)/M << endl;
  Output << radius << endl;
  Output.close();
  
  cout << double(count)/M << endl;

  /*****************************************************/
  
  //initial position
  y[0]=1/sqrt(3);
  y[1]=1/sqrt(3);
  y[2]=1/sqrt(3);
  
  cout << "-------------------------" << endl;
  cout << "Initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  
  //passi e blocchi simulazione
  M=  1E6;
  count=0;
  N=  100;
  L=  M/N;
  
  //grandezze simulazione
  radius=3.8;
  Psi_2_1_0 * psi2 = new Psi_2_1_0(a0);
  
  //equilibrazione della posizione iniziale
  
  Equilibrate_unif (y, psi2, radius, 1000, rnd);
  cout << "Equilibrated initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  cout << "-------------------------" << endl;
  
  //calcolo dei valori della simulazione
  
  Output.open("xyzpsi210.out");
  
  for(int i=0; i<N; i++){
    
    double sum=0;
    
    for(int j=0; j<L; j++){
      for(int k=0; k<y.size(); k++){
        x[k]=y[k];
      }
      
      T_unif(x, radius, rnd);
      double alpha=Eval_Alpha(y, x, psi2);
      Accept(alpha, y, x, count, rnd);
      sum+=g->Eval(sqrt(R2(y)));
      if(j%100==0)Output<< y[0] << "\t" << y[1] << "\t" << y[2] << endl;
      
    }
    
    AV[i] =sum/L;
    AV2[i]=pow(sum/L,2);
    
  }
  
  Output.close();
  Eval_Print("psi210.out", AV, AV2, N);
  
  Output.open("psi210_acceptance.out");
  Output << double(count)/M << endl;
  Output << radius << endl;
  Output.close();
  
  cout << double(count)/M << endl;
  
  /*****************************************************/
  //T-GAUSS
  /*****************************************************/
  
  
  //initial position
  y[0]=1/sqrt(3);
  y[1]=1/sqrt(3);
  y[2]=1/sqrt(3);
  
  cout << "-------------------------" << endl;
  cout << "Initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  
  //passi e blocchi simulazione
  M=  1E6;
  count=0;
  N=  100;
  L=  M/N;
  
  //grandezze simulazione
  double sigma=1.5/2.;
  //equilibrazione della posizione iniziale
  
  Equilibrate_gauss (y, psi1, sigma, 100, rnd);
  cout << "Equilibrated initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  cout << "-------------------------" << endl;
  
  //calcolo dei valori della simulazione
  
  Output.open("xyzpsi100_gauss.out");
  
  for(int i=0; i<N; i++){
    
    double sum=0;
    
    for(int j=0; j<L; j++){
      for(int k=0; k<y.size(); k++){
        x[k]=y[k];
      }
      
      T_gauss(x, sigma, rnd);
      double alpha=Eval_Alpha(y, x, psi1);
      Accept(alpha, y, x, count, rnd);
      sum+=g->Eval(sqrt(R2(y)));
      if(j%100==0)Output<< y[0] << "\t" << y[1] << "\t" << y[2] << endl;
      
    }
    
    AV[i] =sum/L;
    AV2[i]=pow(sum/L,2);
    
  }
  
  Output.close();
  Eval_Print("psi100_gauss.out", AV, AV2, N);
  
  Output.open("psi100_gauss_acceptance.out");
  Output << double(count)/M << endl;
  Output << sigma << endl;
  Output.close();
  
  cout << double(count)/M << endl;
  
  /*****************************************************/
  
  //initial position
  y[0]=1/sqrt(3);
  y[1]=1/sqrt(3);
  y[2]=1/sqrt(3);
  
  cout << "-------------------------" << endl;
  cout << "Initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  
  //passi e blocchi simulazione
  M=  1E6;
  count=0;
  N=  100;
  L=  M/N;
  
  //grandezze simulazione
  sigma=3.8/2.;
  //equilibrazione della posizione iniziale
  
  Equilibrate_gauss (y, psi2, sigma, 100, rnd);
  cout << "Equilibrated initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  cout << "-------------------------" << endl;
  
  //calcolo dei valori della simulazione
  
  Output.open("xyzpsi210_gauss.out");
  
  for(int i=0; i<N; i++){
    
    double sum=0;
    
    for(int j=0; j<L; j++){
      for(int k=0; k<y.size(); k++){
        x[k]=y[k];
      }
      
      T_gauss(x, sigma, rnd);
      double alpha=Eval_Alpha(y, x, psi2);
      Accept(alpha, y, x, count, rnd);
      sum+=g->Eval(sqrt(R2(y)));
      if(j%100==0)Output<< y[0] << "\t" << y[1] << "\t" << y[2] << endl;
      
    }
    
    AV[i] =sum/L;
    AV2[i]=pow(sum/L,2);
    
  }
  
  Output.close(); 
  Eval_Print("psi210_gauss.out", AV, AV2, N);
  
  Output.open("psi210_gauss_acceptance.out");
  Output << double(count)/M << endl;
  Output << sigma << endl;
  Output.close();
  
  cout << double(count)/M << endl;
  
  
  /*****************************************************/
  //T-UNIF - ESPERIMENTO VALORI ALTI E BASSI ACCETTANZA
  /*****************************************************/
  //Accettanza bassa
  
  //initial position
  y[0]=1/sqrt(3);
  y[1]=1/sqrt(3);
  y[2]=1/sqrt(3);
  
  cout << "-------------------------" << endl;
  cout << "Initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  
  //passi e blocchi simulazione
  M=  1E6;
  count=0;
  N=  100;
  L=  M/N;
  
  //grandezze simulazione
  radius=5;
  
  //equilibrazione della posizione iniziale
  
  Equilibrate_unif (y, psi1, radius, 1000, rnd);
  cout << "Equilibrated initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  cout << "-------------------------" << endl;
  
  //calcolo dei valori della simulazione
  
  for(int i=0; i<N; i++){
    
    double sum=0;
    
    for(int j=0; j<L; j++){
      for(int k=0; k<y.size(); k++){
        x[k]=y[k];
      }
      
      T_unif(x, radius, rnd);
      double alpha=Eval_Alpha(y, x, psi1);
      Accept(alpha, y, x, count, rnd);
      sum+=g->Eval(sqrt(R2(y)));
      
    }
    
    AV[i] =sum/L;
    AV2[i]=pow(sum/L,2);
    
  }
  Eval_Print("psi100_low.out", AV, AV2, N);
  
  Output.open("psi100_low_acceptance.out");
  Output << double(count)/M << endl;
  Output << radius << endl;
  Output.close();
  
  cout << double(count)/M << endl;
  
  /*****************************************************/
  
  //initial position
  y[0]=1/sqrt(3);
  y[1]=1/sqrt(3);
  y[2]=1/sqrt(3);
  
  cout << "-------------------------" << endl;
  cout << "Initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  
  //passi e blocchi simulazione
  M=  1E6;
  count=0;
  N=  100;
  L=  M/N;
  
  //grandezze simulazione
  radius=0.2;
  
  //equilibrazione della posizione iniziale
  
  Equilibrate_unif (y, psi1, radius, 1000, rnd);
  cout << "Equilibrated initial position: " << y[0] << "\t" << y[1] << "\t" << y[2] << "\t" << endl;
  cout << "-------------------------" << endl;
  
  //calcolo dei valori della simulazione
  
  for(int i=0; i<N; i++){
    
    double sum=0;
    
    for(int j=0; j<L; j++){
      for(int k=0; k<y.size(); k++){
        x[k]=y[k];
      }
      
      T_unif(x, radius, rnd);
      double alpha=Eval_Alpha(y, x, psi1);
      Accept(alpha, y, x, count, rnd);
      sum+=g->Eval(sqrt(R2(y)));
      
    }
    
    AV[i] =sum/L;
    AV2[i]=pow(sum/L,2);
    
  }
  
  Eval_Print("psi100_high.out", AV, AV2, N);
  
  Output.open("psi100_high_acceptance.out");
  Output << double(count)/M << endl;
  Output << radius << endl;
  Output.close();
  
  cout << double(count)/M << endl;
  
  
  
  return 0;
  }

