#include "random.h"
#include "funzioneBase.h"
#include "variational.h"

using namespace std;


int main (int argc, char *argv[]){

  PRNGen G; //setting generatore
  Random * rnd=G.SetGen("Primes", "seed.in");
  
  double mu=1;
  double sigma=1;
  double rpsi=2.5;
  double rus=0.1;
  double x0=1;
  double ti=1;
  double tf=0.0001;
  
  Psi * psi;
  ELoc * eloc;
  Boltzmann * boltz;
  
  Variational var(psi, eloc, boltz, rnd, mu, sigma, rpsi, rus, x0, ti, tf);
  
  //************************************//
 
  for(int i=0; i<100000; i++) var.MetropolisPSI();
  var.AVE();
  while(var.GetTemp()>tf){
  var.SetValues();
  for(int i=0; i<100000; i++) var.MetropolisPSI();
  var.AVE();
  var.MetropolisSA();
  }
  var.Print();
  for(int i=0; i<100000; i++) var.MetropolisPSI();
  var.PrintSampled("prob.out");
  var.BlockAVE(2000, 100000);
  var.PrintValErrEnergy("energy.final.out");

  return 0;
  }

