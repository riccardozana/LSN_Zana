#include "random.h"
#include "GA.h"

using namespace std;


int main (int argc, char *argv[]){

  PRNGen G; //setting generatore
  Random * rnd=G.SetGen("Primes", "seed.in");
  
  /*****************************************************/
  
  int dim=32;
  int ind=2;
  int step=100000;
  int stepTemp=1000;
  double ti=10;
  double tf=0.0001;
  Population pop(ind, rnd, dim);
  pop.SetRandomPositionSquare();
  pop.Creator(100000);
  pop.Measure();
  pop.InitTemp (ti);
  while(pop.GetTemp()>tf){
    pop.SetTemp(ti, tf);
    pop.Mutation();
    pop.Measure();
    pop.PrintL("Square_L_SA.out");
    pop.Accept();
    pop.Measure();
    pop.PrintL("Square_L_SA.out");
    pop.PrintLBest("Square_Best_SA.out");
  }
  pop.PrintBestPath("Square_BestPath_SA.out");
  
  return 0;
}
  
  

