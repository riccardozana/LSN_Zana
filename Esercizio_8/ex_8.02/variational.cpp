#include "variational.h"


//***********************************************************//
//PRIVATE
//***********************************************************//
double Variational :: TUnif () {
  double dx;
  dx=m_rand->Rannyu(-m_rpsi, m_rpsi);
  return m_x[m_x.size()-1] + dx;
}

void Variational :: TUnifMuSig () {
  double dmu, dsigma;
  
  do{
    dmu=m_rand->Rannyu(-m_rus, m_rus);
    dsigma=m_rand->Rannyu(-m_rus, m_rus);
  }while((dmu*dmu+dsigma*dsigma)>(m_rus*m_rus));
  
  m_mu.push_back(abs(dmu+m_mu[m_mu.size()-1]));
  m_sigma.push_back(abs(dsigma+m_sigma[m_sigma.size()-1]));
}

double Variational :: Error(std::vector<double> AV, std::vector<double> AV2, int n){
  if(n==0){
    return 0;
  }else{
    return sqrt((AV2[n]-pow(AV[n],2))/double(n));
  }
}

//***********************************************************//
//UTILITIES
//***********************************************************//

Variational :: Variational (FunzioneBase2 * psi, FunzioneBase2 * ELoc, FunzioneBase2 * boltz, Random * rnd, double mu, double sigma, double rpsi, double rus, double start, double ti, double tf){
  m_psi=psi;
  m_ELoc=ELoc;
  m_boltzmann=boltz;
  m_rand=rnd;
  m_mu.push_back(mu);
  m_sigma.push_back(sigma);
  m_rpsi=rpsi;
  m_rus=rus;
  m_x0=start;
  m_temp.push_back(ti);
  m_tempf=tf;
  
  m_x.push_back(m_x0);
  m_accepted=0;
  m_attempted=0;
}


void Variational :: Print(){
  ofstream Ene, Temp, Mu, Sigma;
  Ene.open("ene.out");
  Temp.open("temp.out");
  Mu.open("mu.out");
  Sigma.open("sigma.out");
  
  for(int i=0; i<m_E.size(); i++){
    Ene << m_E[i] << endl;
    Temp << m_temp[i] << endl;
    Mu << m_mu[i] << endl;
    Sigma << m_sigma[i] << endl;
  }
  
}

void Variational :: PrintValErrEnergy(string filename){
  ofstream Output;
  Output.open(filename);
  cout << "Acceptance : " << m_accepted/m_attempted << endl;
  
  for(int i=0; i<m_Efinal.size(); i++){
    Output << i << "\t" << m_Efinal[i] << "\t" << m_err[i] << endl;
  }
}

void Variational :: PrintSampled(string filename){
  ofstream Output;
  Output.open(filename);
  for (int i=0; i<m_x.size(); i++) Output << m_x[i] << endl;
}

//***********************************************************//
//CORE ALGORITMHS
//***********************************************************//


void Variational :: MetropolisPSI(){
  
  
  m_psi =new Psi  (m_mu[m_mu.size()-1], m_sigma[m_sigma.size()-1]);
  m_ELoc=new ELoc (m_mu[m_mu.size()-1], m_sigma[m_sigma.size()-1]);
  m_y=TUnif();
  
  double alpha=fmin(1.,(m_psi->Eval(m_y,0.)/m_psi->Eval(m_x[m_x.size()-1],0.)));
  double r=m_rand->Rannyu();
  if(r<alpha){
    m_x.push_back(m_y);
    m_accepted++;
  }
  else m_x.push_back(m_x[m_x.size()-1]);
  m_attempted++;
}


void Variational :: BlockAVE(int L, int M){
  int    N=M/L;
  vector <double> AV(N);
  vector <double> AV2(N);
  m_ELoc=new ELoc (m_mu[m_mu.size()-1], m_sigma[m_sigma.size()-1]);
  
  for(int i=0; i<N; i++){
    
    double sum=0;
    for(int j=0; j<L; j++){
      sum+=m_ELoc->Eval(m_x[i*L+j], 0);
    }
    
    AV[i] =sum/double(L);
    AV2[i]=pow(sum/L,2);
    
  }
  
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
    m_Efinal.push_back(sum_prog[i]);
    m_err.push_back(err_prog[i]);
  }
}
//***********************************************************//
//SIMULATED ANNEALING
//***********************************************************//


void Variational :: MetropolisSA(){
  
  m_boltzmann =new Boltzmann (m_temp[m_temp.size()-1]);
  
  double alpha=fmin(1.,(m_boltzmann->Eval(m_E[m_E.size()-1],m_E[m_E.size()-2])));
  double r=m_rand->Rannyu();
  if(r<alpha){
  }
  else {
    m_E[m_E.size()-1]=m_E[m_E.size()-2];
    m_mu[m_mu.size()-1]=m_mu[m_mu.size()-2];
    m_sigma[m_sigma.size()-1]=m_sigma[m_sigma.size()-2];
  }
}


void Variational :: AVE(){
  
    double sum=0;
    for(int j=0; j<m_x.size(); j++){
      sum+=m_ELoc->Eval(m_x[j], 0);
    }
    m_E.push_back(sum/double(m_x.size()));
    m_x.clear();
    m_x.push_back(m_x0);
}



void Variational ::SetValues(){

  TUnifMuSig();
  double m=0.999;
  m_temp.push_back(m_temp[m_temp.size()-1]*m);
}

