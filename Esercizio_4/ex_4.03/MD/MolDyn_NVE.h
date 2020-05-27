/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
bool restart;
int nSimulations;
int nStepfinal;

//functions
void Input(std::string);
void Move(void);
void ConfFinal(std::string);
void ConfXYZ(int);
void Measure(std::string);
double Force(int, int);
double Pbc(double);
void Restart(std::string);
void Move_WO_Vel(void);
void MeasureR(std::string, int i);
double Error(std::vector<double> , std::vector<double> , int );
void Eval_Print(std::string, std::vector<double>, std::vector<double>, int );
void Eval_Print_SI(std::string, std::vector<double>, std::vector<double>, int, bool);
void MeasureRB(std::string, int, std::vector<std::vector<double> > &);
void FinalRestart(std::string);


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
