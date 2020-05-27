#include "funzioneBase.h"

//class Cos
Cos :: Cos(){
  m_a=0;
  m_b=0;
  m_c=0;
}

Cos :: Cos(double a, double b, double c){
  m_a=a;
  m_b=b;
  m_c=c;
}

Cos :: ~Cos(){}

double Cos :: Eval(double x) const {
  return m_a*cos(m_b*x+m_c);
}

/*****************************************************/

Potenza :: Potenza(FunzioneBase * f, double p){
  m_f=f;
  m_p=p;
}

Potenza :: ~Potenza(){}

double Potenza :: Eval(double x) const {
  return pow(m_f->Eval(x), m_p);
}

/*****************************************************/

TaylorCos :: TaylorCos(){}

TaylorCos :: ~TaylorCos(){}

double TaylorCos :: Eval(double x) const {
  return 1/(1.00452485553)*(M_PI/2-pow(M_PI,3)*pow(x,2)/16+pow(M_PI,5)*pow(x,4)/768); //prefattore di norm
}

/*****************************************************/

Ratio :: Ratio(FunzioneBase * f, FunzioneBase * g){
  m_fu=f;
  m_fd=g;
}

Ratio :: ~Ratio(){}

double Ratio :: Eval(double x) const {
  return m_fu->Eval(x)/m_fd->Eval(x);
}

/*****************************************************/

Product1_2 :: Product1_2(FunzioneBase * f, FunzioneBase2 * g){
  m_f1=f;
  m_f2=g;
}

Product1_2 :: ~Product1_2(){}

double Product1_2 :: Eval(double x, double y) const {
  return m_f1->Eval(x)*m_f2->Eval(x, y);
}

/*****************************************************/

Retta :: Retta(double m, double q){
  m_m=m;
  m_q=q;
}

Retta :: ~Retta(){}

double Retta :: Eval(double x) const {
  return m_m*x+m_q;
}

/*****************************************************/

Psi_1_0_0 :: Psi_1_0_0(double a0){
  m_a0=a0;
}

Psi_1_0_0 :: ~Psi_1_0_0(){}

double Psi_1_0_0 :: Eval(double r, double theta) const {
  double amp= exp(-r/m_a0)/(pow(m_a0, 3./2.)*sqrt(M_PI));
  return pow(amp, 2);
}

/*****************************************************/

Psi_2_1_0 :: Psi_2_1_0(double a0){
  m_a0=a0;
}

Psi_2_1_0 :: ~Psi_2_1_0(){}

double Psi_2_1_0 :: Eval(double r, double theta) const {
  double amp= (pow(m_a0, -5./2.)/8.)*(sqrt(2./M_PI))*r*exp(-r/(2.*m_a0))*cos(theta);
  return pow(amp, 2);
}

/*****************************************************/


