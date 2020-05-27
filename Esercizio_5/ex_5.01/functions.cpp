#include "functions.h"


//Error function: statistical uncertainty estimation
double Error(std::vector<double> AV, std::vector<double> AV2, int n){
   if(n==0){
      return 0;
   }else{
      return sqrt((AV2[n]-pow(AV[n],2))/double(n));
   }
}

//valuta incertezza con metodo a blocchi e stampa valori
void Eval_Print(std::string namefile, std::vector<double> AV, std::vector<double> AV2, int N){
  
  vector<int> x(N);
  iota(x.begin(), x.end(), 1);
  vector<double> sum_prog(N, 0.0);
  vector<double> su2_prog(N, 0.0);
  vector<double> err_prog(N, 0.0);
  
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
  Output.open(namefile);
  if (Output.is_open()){
    for (int i=0; i<N; i++){
      Output << x[i] << "\t" << sum_prog[i]  << "\t" << err_prog[i] << endl;
      }
  } else cerr << "PROBLEM: Unable to open -> " << namefile << endl;
      Output.close();
}

//valuta incertezza con metodo a blocchi e salva valori

void Eval_Store(std::vector<double> &v, std::vector<double> AV, std::vector<double> AV2, int N, int L){
  
  vector<int> x(N);
  iota(x.begin(), x.end(), 1);
  vector<double> sum_prog(N, 0.0);
  vector<double> su2_prog(N, 0.0);
  vector<double> err_prog(N, 0.0);
  
  for (int i=0; i<N; i++){
    for(int j=0; j<i+1; j++){
      sum_prog[i]+= AV[j];
      su2_prog[i]+= AV2[j];
    }
    sum_prog[i]/=(i+1);
    su2_prog[i]/=(i+1);
    err_prog[i]=Error(sum_prog, su2_prog, i );
  }
  
  
  v.push_back(sqrt(sum_prog[N-1])); //carico la radice di r^2
  v.push_back(err_prog[N-1]/(2*sqrt(sum_prog[N-1]))); //propago l'incertezza di r^2 a sqrt(r^2)
}


// calcola raggio al quadrato
double R2(std::vector<double> v){
  double sum=0;
  for(int i=0; i<v.size(); i++){
    sum+=pow(v[i], 2);
  }
  return sum;
}

//Estrae uniformemente xyz sfera centrata in x0
void T_unif (vector<double> & xold, double radius, Random* rnd){
  double x, y, z;
  
  do{
    x=rnd->Rannyu(-radius, radius);
    y=rnd->Rannyu(-radius, radius);
    z=rnd->Rannyu(-radius, radius);
  }while((x*x+y*y+z*z)>(radius*radius));
  
  xold[0]+=x;
  xold[1]+=y;
  xold[2]+=z;
  
}
//Estrae gaussianamente xyz centrato in x0
void T_gauss (vector<double> & xold, double sigma, Random* rnd){
  double x, y, z;
 //do{
   x=rnd->Gauss(0., sigma);
   y=rnd->Gauss(0., sigma);
   z=rnd->Gauss(0., sigma);
 //}while((x*x+y*y+z*z)>(radius*radius));
  
  xold[0]+=x;
  xold[1]+=y;
  xold[2]+=z;
  
}

// valuta alpha per l'algoritmo di metropolis
double Eval_Alpha(vector<double> xold, vector<double> x, FunzioneBase2* p){
  double ry=sqrt(xold[0]*xold[0]+xold[1]*xold[1]+xold[2]*xold[2]);
  double rx=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  double theta_y=acos(xold[2]/ry);
  double theta_x=acos(x[2]/rx);


  double alpha=fmin(1., (p->Eval(rx, theta_x)/p->Eval(ry, theta_y)));
  
  return alpha;
  
}

//decide il passo n+1 dato alpha in metropolis
void Accept (double alpha, vector<double> &xold, vector<double> x, int & count, Random * rnd){
  double r=rnd->Rannyu(0., 1.);
  if(r<=alpha){
    for(int i=0; i<3; i++){
      xold[i]=x[i];
    }
    count++;
  }
}

//mostra equilibrazione del sistema prima di effettuare misurazioni

void Equilibrate_unif_print (vector<double> &xold, FunzioneBase2 * f, double radius, int N, Random * rnd){
  
  ofstream Output;
  Output.open("far_from_origin.out");
  
  vector<double> x(3, 0.0);
  int count=0;
  
  for(int j=0; j<N; j++){
    
    for(int k=0; k<xold.size(); k++){
      x[k]=xold[k];
    }
    
    T_unif(x, radius, rnd);
    double alpha=Eval_Alpha(xold, x, f);
    Accept(alpha, xold, x, count, rnd);
    Output<< xold[0] << "\t" << xold[1] << "\t" << xold[2] << endl;
    
  }
  
  Output.close();
  
}

//equilibra il sistema prima di effettuare misurazioni

void Equilibrate_unif (vector<double> &xold, FunzioneBase2 * f, double radius, int N, Random * rnd){
  
  vector<double> x(3, 0.0);
  int count=0;
  
  for(int j=0; j<N; j++){
    
    for(int k=0; k<xold.size(); k++){
      x[k]=xold[k];
    }
    
    T_unif(x, radius, rnd);
    double alpha=Eval_Alpha(xold, x, f);
    Accept(alpha, xold, x, count, rnd);
  }
}

//equilibra il sistema prima di effettuare misurazioni gauss
void Equilibrate_gauss (vector<double> &xold, FunzioneBase2 * f, double sigma, int N, Random * rnd){
  
  vector<double> x(3, 0.0);
  int count=0;
  
  for(int j=0; j<N; j++){
    
    for(int k=0; k<xold.size(); k++){
      x[k]=xold[k];
    }
    
    T_gauss(x, sigma, rnd);
    double alpha=Eval_Alpha(xold, x, f);
    Accept(alpha, xold, x, count, rnd);
  }
}


