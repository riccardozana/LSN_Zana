#include "GA.h"

//************************************************************//
// INDIVIDUAL //
//************************************************************//


int Individual :: Pbc(int i) // escludo in automatico lo 0 dalle operazioni
{
  if(i >  m_vec.size()) i = i - m_vec.size()+1; //qua ho tolto la conv int
  if(i == m_vec.size()) i = i - m_vec.size()+1;
  if(i == 0)            i = i + m_vec.size()-1;
  if(i < 0)             i = i + m_vec.size()-1;
  return i;
}

Individual :: Individual(int dim, Random * rnd){
  for(int i=0; i<dim; i++) m_vec.push_back(i);
  m_rand=rnd;
}

Individual :: ~Individual(){}


Individual& Individual :: operator=( Individual& other) // Overload
{
  for(int i=0; i<m_vec.size(); i++){
    m_vec[i]=other.GetVector()[i];
  }
  m_L = other.m_L;
  m_rand = other.m_rand;
  m_genbirth = other.m_genbirth;
  
  return *this;
}

//++++++++++++++++++++++++++++++++++++++++++++++//
//UTILITIES - INDIVIDUAL
//++++++++++++++++++++++++++++++++++++++++++++++//

void Individual :: Check (){
  vector<int> sorted;
  for(int i=0; i<m_vec.size(); i++) sorted.push_back(m_vec[i]);
  sort(sorted.begin(), sorted.end());
  for(int i=0; i< m_vec.size()-1; i++){
    if(sorted[i]==sorted[i+1]){
      cout << "Il vettore ha due valori uguali - CHECK" << endl;
      abort();
    }
    if(sorted[m_vec.size()-1]>(m_vec.size()-1)){
      cout << "Il vettore ha valori maggiori della dim - CHECK" << endl;
    }
  }
}


void Individual :: Print(string filename){
  ofstream Output;
  Output.open(filename, ios::app);
  for(int i=0; i<m_vec.size(); i++){
    Output << m_vec[i] << endl;
  }
  Output << endl << endl;
}

//++++++++++++++++++++++++++++++++++++++++++++++//
//MUTATIONS - INDIVIDUAL
//++++++++++++++++++++++++++++++++++++++++++++++//

void Individual :: PairPermutation(int i){
  if(i>=m_vec.size()) {
    cout << "Error i>=N - PAIR" << endl;
    abort();
  }
  
  if(i==0) {
    cout << "Error i=0 PAIR" << endl;
    abort();
  }
  
  int storage=m_vec[i];
  m_vec[i]=m_vec[Pbc(i+1)];
  m_vec[Pbc(i+1)]=storage;
}

void Individual :: Shift(int i, int m, int n){
  if(m>=m_vec.size()-1) {
    cout << "Error m>=N-1 SHIFT" << endl;
    abort();
  }
  if(i>=m_vec.size() or i==0) {
    cout << "Error i>=N or i=0 SHIFT" << endl;
    abort();
  }
  if(n>=m_vec.size()-1) {
    cout << "Error n>=N-1 SHIFT" << endl;
    abort();
  }
  
  if(n+m>=m_vec.size()-1) {
    cout << "Error n+m>=N  SHIFT" << endl;
    abort();
  }
  
  int a[n+m];
  int b[n+m];
  
  for(int j=0; j<m; j++) a[j]=m_vec[Pbc(i+j)];
  for(int j=0; j<n; j++) b[j]=m_vec[Pbc(Pbc(i+m)+j)];
  for(int j=0; j<n; j++) m_vec[Pbc(i+j)]=b[j];
  for(int j=0; j<m; j++) m_vec[Pbc(Pbc(i+n)+j)]=a[j];
  
}

void Individual :: BlockPermutation(int i, int m){
  if(m>=m_vec.size()/2) {
    cout << "Error m>=N/2 BLOCK PERM" << endl;
    abort();
  }
  if(i>=m_vec.size() or i==0) {
    cout << "Error i>=N or i=0 BLOCK PERM" << endl;
    abort();
  }
  
  int a[m];
  int b[m];
  
  for(int j=0; j<m; j++){
    a[j]=m_vec[Pbc(i+j)];
    b[j]=m_vec[Pbc(Pbc(i+m)+j)];
    m_vec[Pbc(Pbc(i+m)+j)]=a[j];
    m_vec[Pbc(i+j)]=b[j];
  }
}

void Individual :: Inversion(int i, int m){
  if(m>=m_vec.size()-1) {
    cout << "Error m>=N-1 INDIV" << endl;
    abort();
  }
  if(i>=m_vec.size() or i==0) {
    cout << "Error i>=N or i=0 INDIV" << endl;
    abort();
  }
  
  vector<int> a;
  for(int j=0; j<m; j++) a.push_back(m_vec[Pbc(i+j)]);
  std::reverse(a.begin(), a.end());
  for(int j=0; j<m; j++) m_vec[Pbc(i+j)]=a[j];
}




//************************************************************//
// POPULATION //
//************************************************************//

Population :: Population(int n, Random * rnd, int dim){
 m_indiv=n;
 m_rand=rnd;
 m_dim=dim;
 m_minL=99999;
 m_aveL=99999;
 for(int i=0; i<m_indiv; i++){
   Individual ind(m_dim, m_rand);
   m_Matr.push_back(ind);
 }
 
 for(int i=0; i<m_indiv; i++){
   m_NextGen.push_back(m_Matr[i]);
 }
}

Population :: ~Population() {}

//++++++++++++++++++++++++++++++++++++++++++++++//
//UTILITIES - POPULATION
//++++++++++++++++++++++++++++++++++++++++++++++//

void Population :: Print(int i, string filename){
  m_Matr[i].Print(filename);
}

void Population :: PrintL(string filename){
  ofstream Output;
  Output.open(filename);
  
  for(int i=0; i<m_indiv; i++){
    Output << m_Matr[i].GetL() << endl;
  }
}

void Population :: PrintLBest(string filename){
  ofstream Output;
  Output.open(filename, ios::app);
  Output << m_gen << "\t" << m_Matr[0].GetL() << endl;
}

void Population :: PrintLAve(string filename){
  ofstream Output;
  Output.open(filename, ios::app);
  double sum=0;
  for(int i=0; i<m_indiv/2; i++) sum+=m_Matr[i].GetL();
  Output << m_gen << "\t" << sum/(m_indiv/2) << endl;
}

void Population :: PrintBestPath(string filename, int rank){
  ofstream Output;
  Output.open(filename);
  cout << "Il migliore della simulazione : " << rank << " è stato ottenuto alla generazione : " << m_best_gen << endl;
  
  for(int i=0; i<m_dim; i++){
    Output << m_BestPath[i] << "\t" << m_pos[0][m_BestPath[i]] << "\t" << m_pos[1][m_BestPath[i]] << endl;
  }
  Output << m_BestPath.size() << "\t" << m_pos[0][0] << "\t" << m_pos[1][0] << endl;
}

void Population :: GlobalBest(int & rank){
  
  int s[m_dim];
  int Sgen, R1gen, R2gen, R3gen;
  double L1, L2, L3;
  int r1[m_dim];
  int r2[m_dim];
  int r3[m_dim];
  int itag1=1;
  int itag2=2;
  int itag3=3;
  int itag4=4;
  int itag5=5;
  int itag6=6;
  int itag7=7;
  int itag8=8;
  int itag9=9;
  
  
  for(int i=0; i<m_dim; i++){
    s[i]=m_BestPath[i];
  }
  Sgen=m_best_gen;
  
  MPI_Status stat1, stat2, stat3, stat4, stat5, stat6, stat7, stat8, stat9;
  
  if(rank==1){
    MPI_Send(s,m_dim,MPI_INTEGER,0,itag1,MPI_COMM_WORLD);
  }else if(rank==0){
    MPI_Recv(r1,m_dim,MPI_INTEGER,1,itag1, MPI_COMM_WORLD,&stat1);
  }
  
  if(rank==0){
    MPI_Recv(r2,m_dim,MPI_INTEGER,2,itag2, MPI_COMM_WORLD,&stat2);
  }else if(rank==2){
   MPI_Send(s,m_dim,MPI_INTEGER,0,itag2,MPI_COMM_WORLD);
  }
  
  if(rank==0){
    MPI_Recv(r3,m_dim,MPI_INTEGER,3,itag3, MPI_COMM_WORLD,&stat3);
  }else if(rank==3){
    MPI_Send(s,m_dim,MPI_INTEGER,0,itag3,MPI_COMM_WORLD);
  }
  
  if(rank==0){
    MPI_Recv(&R1gen,1,MPI_INTEGER,1,itag4, MPI_COMM_WORLD,&stat4);
    MPI_Recv(&L1,1,MPI_DOUBLE,1,itag5, MPI_COMM_WORLD,&stat5);
  }else if(rank==1){
    MPI_Send(&Sgen,1,MPI_INTEGER,0,itag4,MPI_COMM_WORLD);
    MPI_Send(&m_minL,1,MPI_DOUBLE,0,itag5,MPI_COMM_WORLD);
  }

  if(rank==0){
    MPI_Recv(&R2gen,1,MPI_INTEGER,2,itag6, MPI_COMM_WORLD,&stat6);
    MPI_Recv(&L2,1,MPI_DOUBLE,2,itag7, MPI_COMM_WORLD,&stat7);
  }else if(rank==2){
    MPI_Send(&Sgen,1,MPI_INTEGER,0,itag6,MPI_COMM_WORLD);
    MPI_Send(&m_minL,1,MPI_DOUBLE,0,itag7,MPI_COMM_WORLD);

  }
  
  if(rank==0){
    MPI_Recv(&R3gen,1,MPI_INTEGER,3,itag8, MPI_COMM_WORLD,&stat8);
    MPI_Recv(&L3,1,MPI_DOUBLE,3,itag9, MPI_COMM_WORLD,&stat9);
  }else if(rank==3){
    MPI_Send(&Sgen,1,MPI_INTEGER,0,itag8,MPI_COMM_WORLD);
    MPI_Send(&m_minL,1,MPI_DOUBLE,0,itag9,MPI_COMM_WORLD);
  }
  double MIN;
  int GEN; 
  if(rank==0){
    if(fmin(fmin(fmin(m_minL,L1),L2),L3)==L1){
      for(int i=0; i<m_dim; i++){
        m_BestPath[i]=r1[i];
      }
      MIN=L1;
      GEN=R1gen;
    }
    if(fmin(fmin(fmin(m_minL,L1),L2),L3)==L2){
      for(int i=0; i<m_dim; i++){
        m_BestPath[i]=r2[i];
      }
      MIN=L2;
      GEN=R2gen;
    }
    if(fmin(fmin(fmin(m_minL,L1),L2),L3)==L3){
      for(int i=0; i<m_dim; i++){
        m_BestPath[i]=r3[i];
      }
      MIN=L3;
      GEN=R3gen;
    }
    if(fmin(fmin(fmin(m_minL,L1),L2),L3)==m_minL){
      MIN=m_minL;
      GEN=m_best_gen;
    }
    cout << "LUNGHEZZA MINIMA GLOBALE : " << MIN << " GEN " << GEN << endl;
  }
  
  
}

void Population :: CopyBest(int * v){
  for(int i=0; i<m_dim; i++){
    v[i]=m_Matr[0].GetVector()[i];
  }
}

void Population :: SetBest(int * v){
  for(int i=0; i<m_dim; i++){
    m_Matr[0].GetVector()[i]=v[i];
  }
}




//++++++++++++++++++++++++++++++++++++++++++++++//
//SET and MEASURE - POPULATION
//++++++++++++++++++++++++++++++++++++++++++++++//


void Population :: SetPositionCircle (){
  double r=1;
  double theta=0;
  
  double dtheta=2*M_PI/m_dim;
  
  vector<double> x;
  vector<double> y;
  
  for(int i=0; i<m_dim; i++){
    x.push_back(r*cos(theta+dtheta*i));
    y.push_back(r*sin(theta+dtheta*i));
  }
  m_pos.push_back(x);
  m_pos.push_back(y);
  
  ofstream Output;
  Output.open("position_circle.out");
  for(int i=0; i<2; i++){
    for(int j=0; j<m_dim; j++){
      Output << m_pos[i][j]  << endl;
    }
  }
}


void Population :: SetRandomPositionCircle (){
  
  double r=1;
  vector<double> x;
  vector<double> y;
  
  for(int i=0; i<m_dim; i++){
    double theta=m_rand->Rannyu(0, 2*M_PI);
    x.push_back(r*cos(theta));
    y.push_back(r*sin(theta));
  }
  m_pos.push_back(x);
  m_pos.push_back(y);
  
  ofstream Output;
  Output.open("random_position_circle.out");
  for(int i=0; i<2; i++){
    for(int j=0; j<m_dim; j++){
      Output << m_pos[i][j]  << endl;
    }
  }
}

void Population :: SetRandomPositionSquare (){
  
  double l=1;
  vector<double> x;
  vector<double> y;
  
  for(int i=0; i<m_dim; i++){
    x.push_back(m_rand->Rannyu(-l, l));
    y.push_back(m_rand->Rannyu(-l, l));
  }
  m_pos.push_back(x);
  m_pos.push_back(y);
  
  ofstream Output;
  Output.open("random_position_square.out");
  for(int i=0; i<2; i++){
    for(int j=0; j<m_dim; j++){
      Output << m_pos[i][j]  << endl;
    }
  }
}

void Population :: Creator(int n_iter){
  
  for(int i=0; i<n_iter; i++){
    for(int j=0; j<m_indiv; j++) m_Matr[j].PairPermutation(int(m_rand->Rannyu(1, m_dim)));
  }
  
  for(int j=0; j<m_indiv; j++) m_Matr[j].Check();
  m_gen=0;
}

void Population :: Measure (){
  for(int i=0; i<m_indiv; i++){
    double L=0;
    for(int j=0; j<m_dim-1; j++){
      int k1=m_Matr[i].GetVector()[j];
      int k2=m_Matr[i].GetVector()[j+1];
      double dx=m_pos[0][k1]-m_pos[0][k2];
      double dy=m_pos[1][k1]-m_pos[1][k2];
      L+=sqrt(dx*dx+dy*dy);
    }
    int k1=m_Matr[i].GetVector()[m_dim-1];
    int k2=m_Matr[i].GetVector()[0];
    double dx=m_pos[0][k1]-m_pos[0][k2];
    double dy=m_pos[1][k1]-m_pos[1][k2];
    L+=sqrt(dx*dx+dy*dy);
    
    m_Matr[i].SetL(L);
  }
}

void Population :: Sort(){
  
  int min_idx;
  for (int i = 0; i < m_indiv-1; i++){
    min_idx = i;
    for (int j = i+1; j < m_indiv; j++)
      if (m_Matr[j].GetL() < m_Matr[min_idx].GetL()) min_idx = j;
    
    Individual appo = m_Matr[i];
    m_Matr[i] = m_Matr[min_idx];
    m_Matr[min_idx] = appo;
  }
  if(m_Matr[0].GetL()<m_minL){
    m_BestPath=m_Matr[0].GetVector();
    m_minL=m_Matr[0].GetL();
    m_best_gen=m_gen; 
  }
}




//++++++++++++++++++++++++++++++++++++++++++++++//
//MANIPULATION - POPULATION
//++++++++++++++++++++++++++++++++++++++++++++++//


void Population :: Selection(double p){
  
  for(int i=0; i<m_indiv; i=i+2){
    
    int jp=int(m_indiv*pow(m_rand->Rannyu(),p));
    int jm=jp;
    while(jm==jp) jm=int(m_indiv*pow(m_rand->Rannyu(),p));
   // cout << jp << "---" << jm << endl;
    
    m_NextGen[i]=m_Matr[jp];
    m_NextGen[i+1]=m_Matr[jm];
  }
}

void Population :: Crossover(int cut, double p){
  
  for(int i=0; i<m_indiv; i=i+2){
      m_Matr[i]=m_NextGen[i];
      m_Matr[i+1]=m_NextGen[i+1];
      if(m_rand->Rannyu() < p){
        
        int missing1[m_dim-cut];
        int missing2[m_dim-cut];
        for(int j=0; j<m_dim-cut; j++){
          missing1[j]=m_NextGen[i].GetVector()[j+cut];
          missing2[j]=m_NextGen[i+1].GetVector()[j+cut];
        }
        //cout << "debug" << i << endl;
    
        int count1=cut;
        int count2=cut;
        for(int j=1; j<m_dim; j++){ //il primo è sempre 0
          for(int k=0; k<m_dim-cut; k++){
            if(m_NextGen[i+1].GetVector()[j]==missing1[k]){
              m_Matr[i].SetVector(count1, missing1[k]);
              count1++;
            }
            if(m_NextGen[i].GetVector()[j]==missing2[k]){
              m_Matr[i+1].SetVector(count2, missing2[k]);
              count2++;
            }
          }
        }
        for(int k=0; k<m_indiv; k++) m_Matr[k].Check();
      }
    }
  m_gen=m_gen+1;
}

void Population :: Mutation(double p){

  
    for(int i=0; i<m_indiv; i++){
      if(m_rand->Rannyu()<p){
      int j=int(m_rand->Rannyu(1, m_dim));
     /* cout << "Pair" << " ";
      cout << j << endl;*/
      m_Matr[i].PairPermutation(j);
    }
  }
  
    for(int i=0; i<m_indiv; i++){
      if(m_rand->Rannyu()<p){
      int j=int(m_rand->Rannyu(1, m_dim));
      int m=int(m_rand->Rannyu(1,m_dim-1));
      int n=int(m_rand->Rannyu(1,m_dim-1));
      while(m+n>=(m_dim-1)){
        m=int(m_rand->Rannyu(1,m_dim-1));
        n=int(m_rand->Rannyu(1,m_dim-1));
      }
      /*cout << "Shift" << " ";
      cout << j << "\t" << m << "\t" << n << endl;*/
      m_Matr[i].Shift(j, m, n);
    }
   }
  
   for(int i=0; i<m_indiv; i++) m_Matr[i].Check();
  
  
    for(int i=0; i<m_indiv; i++){
      if(m_rand->Rannyu()<p){
      int j=int(m_rand->Rannyu(1, m_dim));
      int m=int(m_rand->Rannyu(1, m_dim/2));
    /*  cout << "Block" << " ";
      cout << j << "\t" << m << endl;*/
      m_Matr[i].BlockPermutation(j, m);
    }
  }
  
  
  if(m_rand->Rannyu()<p){
    for(int i=0; i<m_indiv; i++){
      int j=int(m_rand->Rannyu(1, m_dim));
      int m=int(m_rand->Rannyu(1, m_dim-1));
     /* cout << "inversion" << " ";
      cout << j << "\t" << m << endl;*/
      m_Matr[i].Inversion(j, m);
    }
  }
  for(int i=0; i<m_indiv; i++) m_Matr[i].Check();
  
}



void Population :: SetRelation(int * v, int & target){
  vector<int> c;
  target=m_rand->Rannyu(1,4);
  for(int i=0; i<4; i++) c.push_back(i);
  c.erase(c.begin()+target);
  c.erase(c.begin());
  v[0]=c[0];
  v[1]=c[1];
}



void Population :: Migration(int & rank, int target, int * v){
  int s[m_dim];
  int r[m_dim];
  int itag= 1;
  int itag2=2;
  int itag3=3;
  int itag4=4;
  
  CopyBest(s);
  
  MPI_Status stat1, stat2, stat3, stat4;
 
  if(rank==0){
    cout << rank << "MIGRA VERSO :"  << target << endl;
    MPI_Send(s,m_dim,MPI_INTEGER,target,itag,MPI_COMM_WORLD);
    MPI_Recv(r,m_dim,MPI_INTEGER,target,itag2, MPI_COMM_WORLD,&stat2);
  }else if(rank==target){
    cout << rank << "MIGRA VERSO :"  << 0 << endl;
    MPI_Recv(r,m_dim,MPI_INTEGER,0,itag, MPI_COMM_WORLD,&stat1);
    MPI_Send(s,m_dim,MPI_INTEGER,0,itag2,MPI_COMM_WORLD);
  }
  
  if(rank==v[0]){
    cout << rank << "MIGRA VERSO :"  << v[1] << endl;
    MPI_Send(s,m_dim,MPI_INTEGER,v[1],itag3,MPI_COMM_WORLD);
    MPI_Recv (r,m_dim,MPI_INTEGER,v[1],itag4, MPI_COMM_WORLD,&stat4);
  }else if(rank==v[1]){
    cout << rank << "MIGRA VERSO :"  << v[0] << endl;
    MPI_Recv(r,m_dim,MPI_INTEGER,v[0],itag3, MPI_COMM_WORLD,&stat3);
    MPI_Send(s,m_dim,MPI_INTEGER,v[0],itag4,MPI_COMM_WORLD);
    cout << "---------------" << endl; 
  }
  
  SetBest(r);
}

