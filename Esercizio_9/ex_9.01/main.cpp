#include "random.h"
#include "GA.h"

using namespace std;


int main (int argc, char *argv[]){

  PRNGen G; //setting generatore
  Random * rnd=G.SetGen("Primes", "seed.in");
  
  /*****************************************************/
  
  int dim=32;
  int ind=100;
  Population pop(ind, rnd, dim);
  pop.SetRandomPositionSquare();
  pop.Creator(100000);
  pop.Measure();
  pop.Sort();
  pop.PrintL("Square_InitialL.out");
  
  for(int i=0; i<50000; i++){
    pop.Selection(2);
    pop.Crossover(16, 0.60);
    pop.Mutation(0.07);
    pop.Measure();
    pop.Sort();
    pop.PrintLBest("Square_LBest.out");
    pop.PrintLAve ("Square_LAve.out");
  }
  pop.PrintBestPath("Square_BestPath.out");
  
  return 0;
}
  
  

