#ifndef __variational_h__
#define __variational_h__

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "funzionebase.h"
#include "random.h"


class Variational {
  
private:
  vector <double> m_sigma;
  vector <double> m_mu;
  vector <double> m_Efinal;
  vector <double> m_E;
  vector <double> m_err;
  vector <double> m_x;
  vector <double> m_temp;
  double m_y, m_x0;
  double m_E1;
  double m_accepted, m_attempted, m_rpsi, m_rus, m_tempf;
  Random * m_rand;
  FunzioneBase2 * m_psi;
  FunzioneBase2 * m_ELoc;
  FunzioneBase2 * m_boltzmann;
  
  double TUnif();
  void TUnifMuSig ();
  double Error(std::vector<double> AV, std::vector<double> AV2, int n);


public:

  Variational (FunzioneBase2 * psi, FunzioneBase2 * ELoc, FunzioneBase2 * boltz, Random * rnd, double mu, double sigma, double rpsi, double rus, double start, double ti, double tf);
  ~Variational(){};
  double GetTemp(){return m_temp[m_temp.size()-1];}
  void Print(); 
  void PrintValErrEnergy(string);
  void PrintSampled(string); 

  void MetropolisPSI();
  void BlockAVE(int L, int M);
  
  void MetropolisSA();
  void AVE();
  void SetValues();
  
  
  
};



#endif
