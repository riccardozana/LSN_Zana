 /****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
/*****************************************************************/
  
  Input(); //Inizialization
  PrintConfig("config.0"); //stampo configurazione iniziale per confrontarla
  
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);     //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
    }
    if (restart==0)
    {
      Averages(iblk, 0); //Print results for current block
    }
  }
  if (restart==0) ConfFinal("config/config0.final"); //Write final configuration
  else
  {
    ConfFinal("config/config0.old");
  
    Input("config/config0.old");   //Inizialization
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
        Move();
        Measure();
        Accumulate(); //Update block averages
      }
        Averages(iblk, 0); //Print results for current block
    }
    ConfFinal("config/config0.final");
  }
/*****************************************************************/
  
  
  double dt=(tempMax-temp)/double(tempstep);
  
  for (int i=1; i<tempstep+1; i++){
    Input(dt);
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
        Move();
        Measure();
        Accumulate(); //Update block averages
      }
      if (restart==0)
      {
        Averages(iblk, i); //Print results for current block
      }
    }
    if (restart==0) ConfFinal("config/config" + to_string(i) + ".final"); //Write final configuration
    else
    {
      ConfFinal("config/config" + to_string(i) + ".old");
    
      Input("config/config" + to_string(i) + ".old"); //Inizialization
      for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
      {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
          Move();
          Measure();
          Accumulate(); //Update block averages
        }
          Averages(iblk, i); //Print results for current block
      }
      ConfFinal("config/config" + to_string(i) + ".final");
    } // end else
  }
return 0;
}

/*****************************************************************/

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;
  
  ReadInput >> tempMax;
  
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Input(std::string filename)
{
  cout << "-- Restart from a previous configuration --" << endl;
  
  ifstream ReadConfig;
  ReadConfig.open(filename);
  
  //initial configuration
  for (int i=0; i<nspin; ++i) ReadConfig >> s[i];
  ReadConfig.close();
  
  accepted=0;
  attempted=0;
  
}

void Input(double dt){
  
  //initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
  temp=temp+dt;
  beta=1./temp;
  accepted=0;
  attempted=0;
  
  cout << endl << "NEW SIMULATION : " << endl;
  cout << "TEMPERATURE : " << temp << endl << endl;
  
}

void Move()
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
      double alpha=fmin(1., exp(-1*beta*Ediff(s[o], o)));
      Accept(s, rnd, o, alpha);
      // INCLUDE YOUR CODE HERE***************************/
      attempted = attempted + 1.0;
    }
    else //Gibbs sampling
    {
      s[o]=Gibbs_sample(s, rnd, o);
      
      // INCLUDE YOUR CODE HERE***************************/
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0, c=0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
    //c += pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]), 2);
    m += s[i];
    
// INCLUDE YOUR CODE HERE
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[ix] = pow(m, 2);
  walker[im] = m;
  
  
// INCLUDE YOUR CODE HERE
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk, int n) //Print results for current block
{
    
  ofstream Ene, Heat, Mag, Chi, OutAcceptance;
  const int wd=12;
  int out_h;
  if(h==0) out_h=0;
  else out_h=1;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
  
  if(iblk==20 and h==0. and metro==1)
  {
    OutAcceptance.open("output/Acceptance" + to_string(metro) + to_string(out_h) + ".out", ios::app);
    OutAcceptance << setw(wd) << accepted/attempted << endl;
    OutAcceptance.close();
  }
    
  Ene.open("output/output" + to_string(n) + "." + to_string(metro) + to_string(out_h) + ".ene", ios::app);
  stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);
  Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
  Ene.close();

  Heat.open("output/output" + to_string(n) + "." + to_string(metro) + to_string(out_h) + ".heat", ios::app);
  stima_c = beta*beta*(blk_av[ic]/blk_norm-pow(blk_av[iu]/blk_norm,2))/(double)nspin; //Heat
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);
  Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  Heat.close();

  Chi.open("output/output" + to_string(n) + "." + to_string(metro) + to_string(out_h) + ".chi", ios::app);
  stima_x = beta *( blk_av[ix]/blk_norm/(double)nspin ); //Chi
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);
  Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  Chi.close();

  Mag.open("output/output" + to_string(n) + "." + to_string(metro) + to_string(out_h) + ".mag", ios::app);
  stima_m = ( blk_av[im]/blk_norm/(double)nspin ); //Chi
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);
  Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
  Mag.close();

// INCLUDE YOUR CODE HERE

  cout << "----------------------------" << endl << endl;
}


void ConfFinal(std::string filename)
{
  ofstream WriteConf;

  cout << "Print final configuration to file " + filename << endl << endl;
  WriteConf.open(filename);
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/***************************************************************/

double Ediff(int sm, int ip)
{
  double diff_ene = 2.*J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)]) +2.*sm*h;
  return diff_ene;
}

void Accept(double * s, Random & rnd, int ip, double alpha){
  double r=rnd.Rannyu();
  if(r<alpha)
  {
    s[ip]=-1*s[ip];
    accepted = accepted + 1.0;
  }
}

double Gibbs_sample(double * s, Random & rnd, int ip){
  double r=rnd.Rannyu();
  if(s[Pbc(ip-1)] + s[Pbc(ip+1)] == 0)
  {
    double p1=1/(1+exp(-2.*beta*h));
    if(r<p1) return 1;
    else      return -1;
  }
  else
  {
    double  p1=1/(1+exp(-2.*beta*J*(s[Pbc(ip-1)] + s[Pbc(ip+1)])-2*h));
    if(r<p1) return 1;
    else     return -1;
  }
}

void PrintConfig(std::string filename){
  ofstream WriteConf;
  
  cout << "Print configuration to file " + filename << endl << endl;
  WriteConf.open(filename);
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

}

  



/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
