#include "random.h"
#include "GA.h"
#include "mpi.h"

using namespace std;


int main (int argc, char *argv[]){

  int dim=32;
  int ind=100;
  int Nmigr=200;
  int Gen=4000;
  
  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  /*****************************************************/
  PRNGen G; //setting generatore
  Random * rnd;
  if(rank==0)rnd=G.SetGen("Primes", "seed.in", 0);
  if(rank==1)rnd=G.SetGen("Primes", "seed.in", 1);
  if(rank==2)rnd=G.SetGen("Primes", "seed.in", 2);
  if(rank==3)rnd=G.SetGen("Primes", "seed.in", 3);
  /*****************************************************/
  Population pop(ind, rnd, dim);
  
  //broadcast della posizione iniziale agli altri
  pop.SetRandomPositionSquare();
  double x_pos[dim];
  double y_pos[dim];
  
  if(rank==0){
    for(int i=0; i<dim; i++){
      x_pos[i]=pop.GetPos(0, i);
      y_pos[i]=pop.GetPos(1, i);
    }
  }
  
  MPI_Bcast(x_pos,32,MPI_REAL8,0, MPI_COMM_WORLD);
  MPI_Bcast(y_pos,32,MPI_REAL8,0, MPI_COMM_WORLD);
  
  if(rank!=0){
    for(int i=0; i<dim; i++){
      pop.SetPos(0, i, x_pos[i]);
      pop.SetPos(1, i, y_pos[i]);
    }
  }
  //*******************************************//
  //evoluzione con migrazione
  
  pop.Creator(100000);
  pop.Measure();
  pop.Sort();
  if(rank==0)pop.PrintL("Square_InitialL.0.out");
  if(rank==1)pop.PrintL("Square_InitialL.1.out");
  if(rank==2)pop.PrintL("Square_InitialL.2.out");
  if(rank==3)pop.PrintL("Square_InitialL.3.out");
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(int i=0; i<Gen; i++){
    pop.Selection(2);
    pop.Crossover(16, 0.60);
    pop.Mutation(0.07);
    pop.Measure();
    pop.Sort();
    if(rank==0)pop.PrintLBest("Square_LBest.0.out");
    if(rank==0)pop.PrintLAve ("Square_LAve.0.out");
    if(rank==1)pop.PrintLBest("Square_LBest.1.out");
    if(rank==1)pop.PrintLAve ("Square_LAve.1.out");
    if(rank==2)pop.PrintLBest("Square_LBest.2.out");
    if(rank==2)pop.PrintLAve ("Square_LAve.2.out");
    if(rank==3)pop.PrintLBest("Square_LBest.3.out");
    if(rank==3)pop.PrintLAve ("Square_LAve.3.out");
    if((i+1)%Nmigr==0){
      int v[2];
      int target;
      if(rank==0) pop.SetRelation(v, target);
      MPI_Bcast(&target,1,MPI_INTEGER,0, MPI_COMM_WORLD);
      MPI_Bcast(v,2,MPI_INTEGER,0, MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      pop.Migration(rank, target, v);
    }
  }
  
  if(rank==0)pop.PrintBestPath("Square_BestPath.0.out", 0);
  if(rank==1)pop.PrintBestPath("Square_BestPath.1.out", 1);
  if(rank==2)pop.PrintBestPath("Square_BestPath.2.out", 2);
  if(rank==3)pop.PrintBestPath("Square_BestPath.3.out", 3);
  
  pop.GlobalBest(rank);
  if(rank==0) pop.PrintBestPath("Square_GlobalBestPath.out", 0123);
  MPI_Finalize();
  return 0;
}
  
  

