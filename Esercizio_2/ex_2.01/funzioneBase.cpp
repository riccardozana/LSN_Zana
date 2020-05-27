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

Retta :: Retta(double m, double q){
  m_m=m;
  m_q=q;
}

Retta :: ~Retta(){}

double Retta :: Eval(double x) const {
  return m_m*x+m_q;
}

/*****************************************************/


