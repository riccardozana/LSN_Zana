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
#include <vector>
#include "Monte_Carlo_NVT.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  //equilibration phase
    for(int istep=1; istep <= eqstep; ++istep)
    {
      Move();
      Measure();
      inst_val.push_back(walker[iv]/(double)npart + vtail);
      inst_val.push_back(rho * temp + (walker[iw] + (double)npart * ptail) / vol);
      
      if(istep%100==0)
      {
        cout << "Equilibration phase : " << istep/100 << " of " << eqstep/100 << endl;
      }
    }
  
  Print("instval." + phase + ".out", inst_val);
  inst_val.clear();
  if(inst==1){
  // Valori instantanei equilibrati per script python
  for(int istep=1; istep <= 500000; ++istep)
  {
    Move();
    Measure();
    inst_val.push_back(walker[iv]/(double)npart + vtail);
    inst_val.push_back(rho * temp + (walker[iw] + (double)npart * ptail) / vol);
    
    if(istep%10000==0)
    {
      cout << "Stampo valori istantanei : " << istep/10000 << " of " << 500000/10000 << endl;;
    }
  }
  
  Print("instval_equilibrated." + phase + ".out", inst_val);
  }
  //evaluation phase
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      if(istep%10 == 0){
     //   ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal();  //Write final configuration
  

  return 0;
}


void Input(void)
{
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

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

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl; 

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;
  
  ReadInput >> eqstep;
  
  ReadInput >> phase;
  
  ReadInput >> inst;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
 
  n_props = 2; //Number of observables

//measurement of g(r)
  igofr = 2;
  nbins = 200;
  nbins_f = 10;
  nbins_s = nbins/nbins_f;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;
  bin_f_size = (box/2.0)/(double)nbins_f;
  bin_s_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i)
  {
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = Pbc( x[i] * box );
    y[i] = Pbc( y[i] * box );
    z[i] = Pbc( z[i] * box );
  }
  ReadConf.close();
  
//Evaluate potential energy and virial of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}


void Move(void)
{
  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;


  for(int i=0; i<npart; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(rnd.Rannyu()*npart);

  //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];

    energy_old = Boltzmann(xold,yold,zold,o);

  //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    energy_new = Boltzmann(xnew,ynew,znew,o);

  //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu())  
    {
    //Update
       x[o] = xnew;
       y[o] = ynew;
       z[o] = znew;
    
       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
  }
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

void Measure()
{
  int bin;
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {

// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
      
     
      for(int k=1; k<nbins+1; k++){
        if(dr<k*bin_size){
          walker[igofr+(k-1)]+=2;
          break;
        }
      }
      
      /*for(int k=1; k<nbins_f+1; k++){
        if(dr<k*bin_f_size){
          for(int s=1; s<nbins_s+1; s++){
            if(dr<((k-1)*bin_f_size+s*bin_s_size)){
              walker[igofr+(k-1)*int(nbins/nbins_f)+(s-1)]+=2.;
             // cout << igofr+(k-1)*int(nbins/nbins_f)+(s-1) << endl;
              break;
            }
          }
          break;
        }
      }*/
      
      
      
//update of the histogram of g(r)

     if(dr < rcut)
     {
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

// contribution to energy and virial
       v += vij;
       w += wij;
     }
    }          
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;
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


void Averages(int iblk) //Print results for current block
{
    
   double r, gdir;
   ofstream Gofr, Gave, Epot, Pres;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Epot.open("output.epot." + phase + ".0",ios::app);
    Pres.open("output.pres." + phase + ".0",ios::app);
    Gofr.open("output.gofr." + phase + ".0",ios::app);
    Gave.open("output.gave." + phase + ".0",ios::app);
    
    stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
    
  //  stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
  stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart )/vol;
  cout << vol << endl;
  cout << blk_norm << endl;
  cout << rho * temp << endl;
  cout << ptail * (double)npart << endl;
  cout << epsilon/pow(sigma,3);
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_press=Error(glob_av[iw],glob_av2[iw],iblk);
  

//Potential energy per particle
    Epot << iblk <<  "\t" << stima_pot*epsilon << "\t" << glob_av[iv]/(double)iblk*epsilon << "\t" << err_pot*epsilon << endl;
//Pressure
   Pres << iblk <<  "\t" << stima_pres*epsilon/pow(sigma,3) << "\t" << glob_av[iw]/(double)iblk*epsilon/pow(sigma,3) << "\t" << err_press*epsilon/pow(sigma,3) << endl;

//g(r)
  

  for(int i=0; i<nbins-1; i++){
    // valori medi in funzione di r
    Gofr << "---------Blocco : " << i;
    Garr[i]=blk_av[igofr+i]/blk_norm/(DeltaV(bin_size*i)*rho*npart);
    Gofr << Garr[i] << setw(wd);
  }
  Gofr << endl << endl;
  
  for(int i=0; i<nbins-1; i++){
    glob_av[igofr+i] += Garr[i];
    glob_av2[igofr+i] += Garr[i]*Garr[i];
    if(iblk==nblk){
      err_gdir=Error(glob_av[igofr+i],glob_av2[igofr+i],iblk);
      Gave << (bin_size*i) << "\t" << glob_av[igofr+i]/(double)iblk << "\t" << err_gdir << endl;
    }
  }
  
  

  cout << "----------------------------" << endl << endl;

  Epot.close();
  Pres.close();
  Gofr.close();
  Gave.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config." + phase + ".final");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/***************************************************************/

void Print(std::string filename, std::vector <double> v) {
  ofstream Output;
  Output.open(filename);
  for(int i=0; i<v.size(); i=i+2){
    Output << v[i] << "\t" << v[i+1] << endl;
  }
  Output.close();
}

double DeltaV(double r){
  return (4.*M_PI)/3*(pow(r + bin_size, 3)-pow(r, 3));
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
