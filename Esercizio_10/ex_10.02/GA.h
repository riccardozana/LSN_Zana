#ifndef __GA_h__
#define __GA_h__

#include <cmath>
#include "random.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include "mpi.h"

class Individual{

public:
  //overload
  Individual& operator=(Individual&);
  //costructors
  Individual(int, Random *);
  // destructor
  ~Individual();
  
  void Check();
  void Print(string);
  void SetL(double L){m_L=L;}
  
  vector<int> GetVector(){return m_vec;}
  void SetVector(int i, int value){m_vec[i]=value;}
  double GetL(){return m_L;}
  
  //Mutations
  void PairPermutation(int);
  void Shift(int , int , int);
  void BlockPermutation(int, int);
  void Inversion(int, int);

private:
  vector<int> m_vec; //vettore con città
  int         m_genbirth; //generazione
  Random *    m_rand;
  double      m_L;
  int         Pbc(int i);

};

class Population{
  
private:
  int                     m_indiv, m_gen, m_dim, m_best_gen;
  double                  m_minL, m_aveL;
  Random *                m_rand;
  vector<Individual>      m_Matr;
  vector<Individual>      m_NextGen;
  vector<vector<double> > m_pos;
  vector<int>             m_BestPath;

  
  
  int Pbc(int i);  //Algorithm for periodic boundary conditions
  
  
public:
  //costructors
  Population(int n, Random * rnd, int dim);
  // destructor
  ~Population();
  
  void   Print(int, string);
  void   PrintL(string);
  void   PrintLBest(string);
  void   PrintLAve(string);
  void   Copy(int i, int j){m_Matr[i]=m_Matr[j];}
  void   PrintBestPath(string, int);
  void   GlobalBest(int &);
  double GetPos(int i, int j){return m_pos[i][j];}
  void   SetPos(int i, int j, double val){m_pos[i][j]=val;}
  void   CopyBest(int *);
  void   SetBest(int *);
  
  
  void SetPositionCircle();
  void SetRandomPositionCircle();
  void SetRandomPositionSquare();
  void Creator(int);
  void Measure();
  void Sort();
 // void Selection(double);
  
  void Selection(double);
  void Crossover(int, double);
  void Mutation(double);
  void SetRelation(int * v, int & target);
  void Migration(int & rank, int target, int * v);
};
#endif


