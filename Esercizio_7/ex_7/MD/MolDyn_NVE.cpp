/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>
#include <ostream>
#include <fstream>
#include <cmath>// rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input();
  
  for(int i=0; i<nSimulations; i++)
  {
    cout << endl << endl << "Restart N. " << i+1 << " of " << nSimulations << endl << endl;
    for(int istep=0; istep<=eqstep; istep++)
    {
      Move();
      if(istep%iprint == 1000) cout << "Equilibration Phase - Number of time-steps: " << istep << endl;
      Measure();
    }
    
    //Riscalo velocità
    //calcolo v(t+dt/2)
    for(int i=0; i<npart; ++i)
    {
      vx[i] = Pbc(x[i] - xold[i])/(delta);
      vy[i] = Pbc(y[i] - yold[i])/(delta);
      vz[i] = Pbc(z[i] - zold[i])/(delta);
    }
    
    //stimo Temp(t+dt/2)
    double t = 0.0;
    for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    stima_temp = (2.0 / 3.0) * t/(double)npart;
    cout << "Temp (t+dt/2): " << stima_temp << endl;
    
    //riscalo le velocità
    double c=temp/stima_temp;
    for(int i=0; i<npart; ++i)
    {
      vx[i] = sqrt(c)*vx[i];
      vy[i] = sqrt(c)*vy[i];
      vz[i] = sqrt(c)*vz[i];
    }
    
    //stimo la nuova configurazione r(t)
    
    for(int i=0; i<npart; ++i)
    {
      xold[i] = Pbc( x[i] - vx[i] * delta );
      yold[i] = Pbc( y[i] - vy[i] * delta );
      zold[i] = Pbc( z[i] - vz[i] * delta );
    }

  }
  
  //final equilibration phase
  temp=stima_temp;
  for(int istep=0; istep<=10000; istep++)
  {
    Move();
    if(istep%iprint == 1000) cout << "Final equilibration phase..." << endl;
  }
  
  //Evaluation phase
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk)
  {
    Reset(iblk);
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0) Measure();
      if(istep%10 == 0) Accumulate();
      if(istep%10 == 0) Print_temp("temp_insta.out");
      
      if(istep%10 == 0){
        //   ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
      nconf += 1;
      }
    }
    Averages(iblk);
  }
  ConfFinal();
  
  return 0;
}

/**********************************************************/

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 2;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> eqstep;
  ReadInput >> phase;
  ReadInput >> iprint;
  ReadInput >> nSimulations;
  
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4;
  n_props = 4; //Number of observables
  
  //measurement of g(r)
  igofr = 5;
  nbins = 200;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v = 0.0, w = 0.0, t=0.0;
  double vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Vir;
  
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  Epot.open(phase + "_out/output_epot.dat",ios::app);
  Ekin.open(phase + "_out/output_ekin.dat",ios::app);
  Temp.open(phase + "_out/output_temp.dat",ios::app);
  Etot.open(phase + "_out/output_etot.dat",ios::app);
  Vir.open(phase + "_out/output_vir.dat",ios::app);

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
      
    for(int k=1; k<nbins+1; k++){
      if(dr<k*bin_size){
        walker[igofr+(k-1)]+=2;
        break;
      }
    }

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
       
//Potential energy
       v += vij;
       w += wij;
     }
    }          
  }

  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  
    walker[ik]=t;
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Vir  << w << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Vir.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
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

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
  ofstream Gofr, Gave, Epot, Pres, Ekin, Etot, Temp;
  const int wd=10;
  
  cout << "Block number " << iblk << endl;
  
  Epot.open(phase + "_out/ave_epot.out",ios::app);
  Pres.open(phase + "_out/ave_press.out",ios::app);
  Ekin.open(phase + "_out/ave_ekin.out",ios::app);
  Etot.open(phase + "_out/ave_etot.out",ios::app);
  Temp.open(phase + "_out/ave_temp.out",ios::app);
  Gofr.open(phase + "_out/gofr.out",ios::app);
  Gave.open(phase + "_out/ave_gave.out",ios::app);
  
  stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
  
  //stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
  stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail*(double)npart) / vol;
  cout << vol << endl;
  cout << blk_norm << endl;
  cout << rho * temp << endl;
  cout << ptail*(double)npart << endl;
  cout << epsilon/pow(sigma,3);
  glob_av[iw] += stima_pres;
  glob_av2[iw] += stima_pres*stima_pres;
  err_press=Error(glob_av[iw],glob_av2[iw],iblk);
  
  stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic Energy
  glob_av[ik] += stima_kin;
  glob_av2[ik] += stima_kin*stima_kin;
  err_kin=Error(glob_av[ik],glob_av2[ik],iblk);
  
  stima_etot = (stima_kin+stima_pot); //Total Energy
  glob_av[ie] += stima_etot;
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot=Error(glob_av[ie],glob_av2[ie],iblk);
  
  stima_temp = (2.0 / 3.0) * blk_av[ik]/blk_norm/(double)npart; //Temp
  glob_av[it] += stima_temp;
  glob_av2[it] += stima_temp*stima_temp;
  err_temp=Error(glob_av[it],glob_av2[it],iblk);
  
  //Potential energy per particle
  Epot << iblk <<  "\t" << stima_pot*epsilon << "\t" << glob_av[iv]/(double)iblk*epsilon << "\t" << err_pot*epsilon << endl;
  //Pressure
 Pres << iblk <<  "\t" << stima_pres*epsilon/pow(sigma,3) << "\t" << glob_av[iw]/(double)iblk*epsilon/pow(sigma,3) << "\t" << err_press*epsilon/pow(sigma,3) << endl;
  //Kinetic energy per particle
  Ekin << iblk <<  "\t" << stima_kin*epsilon << "\t" << glob_av[ik]/(double)iblk*epsilon << "\t" << err_kin*epsilon << endl;
  //Total energy per particle
  Etot << iblk <<  "\t" << stima_etot << "\t" << glob_av[ie]/(double)iblk << "\t" << err_etot << endl;
  //Temperature
  Temp << iblk <<  "\t" << stima_temp*epsilon/boltzmann << "\t" << glob_av[it]/(double)iblk*epsilon/boltzmann << "\t" << err_temp*epsilon/boltzmann << endl;
  
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


double Error(double sum, double sum2, int iblk)
{
  if( iblk == 1 ) return 0.0;
  else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

double DeltaV(double r){
  return (4.*M_PI)/3*(pow(r + bin_size, 3)-pow(r, 3));
}

void Print_temp(std::string filename){
  ofstream Temp;
  Temp.open(phase + "_out/" + filename, ios::app);
  Temp << stima_temp << endl;
  Temp.close();
}
