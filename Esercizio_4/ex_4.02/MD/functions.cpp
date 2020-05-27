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

//valutazione S(t) diretta
double S_direct(double S_0, double mu, double sgm, double T, Random * rnd){
  double w=rnd->Gauss(0, T);
  return S_0*exp((mu-0.5*pow(sgm,2))*T+(sgm*w));
}

//valutazione S(t) discreta
double S_discret(double S_0, double mu, double sgm, double T, Random * rnd, int N){
  double S_ti=S_0;
  double S_tf;
  double dt=T/N;
  
  for(int i=0; i<N; i++){
    double zi=rnd->Gauss(0, 1);
    S_tf=S_ti*exp((mu-0.5*pow(sgm,2))*dt+sgm*zi*sqrt(dt));
    S_ti=S_tf; 
  }
 
  return S_tf;
}

double C00(double S_T, double r, double T,  double K){
  return exp(-r*T)*fmax(0, S_T-K);
}

double P00(double S_T, double r, double T,  double K){
  return exp(-r*T)*fmax(0, K-S_T);

}

