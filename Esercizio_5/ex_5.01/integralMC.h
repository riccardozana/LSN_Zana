#ifndef __IntegralMC_h__
#define __IntegralMC_h__

#include "random.h"
#include "funzioneBase.h"
#include <cmath>

using namespace std;


class IntegralMC {

private:
  Random * m_myrand;
  
protected:
  
public:
	// constructor
  IntegralMC(Random * rnd) {m_myrand = rnd ;};
  // destructor
	~IntegralMC() {;};
  // methods
  double IntegralAVE(double xmin, double xmax, FunzioneBase * f, int punti) ;
  double RejTec(double xmin, double xmax, double ymax, FunzioneBase * f);
  double RejTec_Hybrid(double xmin, double xmax, FunzioneBase * g,  FunzioneBase * p);
  double IntegralIMP(double xmin, double xmax, double ymax, FunzioneBase * p, FunzioneBase * f, int punti);
  double Integral_IMP_Hybrid(double xmin, double xmax, FunzioneBase * g, FunzioneBase * p, FunzioneBase * f, int punti);
};

#endif //__IntegralMC_h__





