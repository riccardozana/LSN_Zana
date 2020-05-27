#include "integralMC.h"

double IntegralMC::IntegralAVE(double xmin, double xmax, FunzioneBase * f, int punti){
	
  double x=0;
	double sum=0;
  
  for(int i=0; i<punti; i++){
    double r=m_myrand->Rannyu(xmin, xmax);
    sum+=f->Eval(r);}

	return sum*(xmax-xmin)/punti;
}

double IntegralMC::RejTec(double xmin, double xmax, double ymax, FunzioneBase * f){
  
  double x;
  double y;
  do{
    x=m_myrand->Rannyu(xmin, xmax);
    y=m_myrand->Rannyu();
  }while(y>(f->Eval(x)/ymax));
  
  return x;
}

double IntegralMC::RejTec_Hybrid(double xmin, double xmax, FunzioneBase * g,  FunzioneBase * p){
  
  double x;
  double y;
  do{
    x=m_myrand-> G_xHybrid();
    y=m_myrand->Rannyu();
  }while(y>(p->Eval(x)/g->Eval(x)));
  
  return x;
}



double IntegralMC::IntegralIMP(double xmin, double xmax, double ymax, FunzioneBase * p, FunzioneBase * f, int punti){
  double x=0;
  double sum=0;
  
  for(int i=0; i<punti; i++){
    double r=RejTec(xmin, xmax, ymax, p);
    sum+=f->Eval(r);}
  
  return sum*(xmax-xmin)/punti;
  
}

double IntegralMC::Integral_IMP_Hybrid(double xmin, double xmax, FunzioneBase * g, FunzioneBase * p, FunzioneBase * f, int punti){
  double x=0;
  double sum=0;
  
  for(int i=0; i<punti; i++){
    double r=RejTec_Hybrid(xmin, xmax, g, p);
    sum+=f->Eval(r);}
  
  return sum*(xmax-xmin)/punti;
  
}








