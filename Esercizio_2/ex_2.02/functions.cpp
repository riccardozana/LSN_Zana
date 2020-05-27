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
void Eval_Print(std::string namefile, std::vector<double> AV, std::vector<double> AV2, int N, int L){
  
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
      Output << x[i]*L << "\t" << sum_prog[i]  << "\t" << err_prog[i] << endl;
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

