#ifndef __funzioneBase_h__
#define __funzioneBase_h__

#include <cmath>

using namespace std;

class FunzioneBase{

public:
	virtual double Eval(double) const =0;
	
};

class Cos : public FunzioneBase {
  
private:
  double m_a, m_b, m_c;
  
protected:

public:
  // constructors
  Cos();
  Cos(double a, double b, double c);
  // destructor
  ~Cos();
  //  methods
  double GetA() const {return m_a;};
  double GetB() const {return m_b;};
  double GetC() const {return m_c;};
  
  void SetA(double a) {m_a=a;};
  void SetB(double b){m_b=b;};
  void SetC(double c) {m_c=c;};
  
  virtual double Eval(double x) const;

};

/*****************************************************/

class Potenza : public FunzioneBase {
  
private:
  FunzioneBase* m_f;
  double m_p;
  
protected:
  
public:
  // constructors
  Potenza(FunzioneBase * f, double p);
  // destructor
  ~Potenza();
  //  methods
  virtual double Eval(double x) const;
  
};

/*****************************************************/

class TaylorCos : public FunzioneBase {
  
private:
  
protected:
  
public:
  // constructors
  TaylorCos();
  // destructor
  ~TaylorCos();
  //  methods
  virtual double Eval(double x) const;
  
};

/*****************************************************/

class Ratio : public FunzioneBase {
  
private:
  FunzioneBase * m_fu;
  FunzioneBase * m_fd;
  
protected:
  
public:
  // constructors
  Ratio(FunzioneBase *, FunzioneBase *);
  // destructor
  ~Ratio();
  //  methods
  virtual double Eval(double x) const;
  
};

/*****************************************************/

class Retta : public FunzioneBase {
  
private:
  double m_m;
  double m_q;
  
protected:
  
public:
  // constructors
  Retta(double m, double q);
  // destructor
  ~Retta();
  //  methods
  virtual double Eval(double x) const;
  
};

/*****************************************************/
/*****************************************************/
/*****************************************************/

class FunzioneBase2{
  
public:
  virtual double Eval(double, double) const =0;
  
};

/*****************************************************/

class Psi_1_0_0 : public FunzioneBase2 {
  
private:
  double m_a0;
  
  
protected:
  
public:
  // constructors
  Psi_1_0_0(double a0);
  // destructor
  ~Psi_1_0_0();
  //  methods
  virtual double Eval(double r, double theta) const;
  
};

/*****************************************************/

class Psi_2_1_0 : public FunzioneBase2 {
  
private:
  double m_a0;
  
  
protected:
  
public:
  // constructors
  Psi_2_1_0(double a0);
  // destructor
  ~Psi_2_1_0();
  //  methods
  virtual double Eval(double r, double theta) const;
  
};

/*****************************************************/

class Product1_2 : public FunzioneBase2 {
  
private:
  FunzioneBase * m_f1;
  FunzioneBase2 * m_f2;
  
protected:
  
public:
  // constructors
  Product1_2(FunzioneBase *, FunzioneBase2 *);
  // destructor
  ~Product1_2();
  //  methods
  virtual double Eval(double x, double y) const;
  
};


/*****************************************************/






#endif

