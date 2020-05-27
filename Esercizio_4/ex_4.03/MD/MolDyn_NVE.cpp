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
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <numeric>
#include "MolDyn_NVE.h"


using namespace std;

//definisco i valori per la conversione come variabili globali
double sigma=0.34; //nm
double boltzmann =1.38064852E-23; //J/k
double epsilon=boltzmann*120; //J
double m=39.948; //uma

int main(){

  
  //SIMULAZIONE FASE SOLIDA //
  
  Input("input_complete.solid");
  for(int istep=1; istep <= nstep; ++istep){
     Move();
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure("solid_out");
     }
  }
  Restart("solid_out");
  FinalRestart("solid_out");
  
  //SIMULAZIONE FASE LIQUIDA
  
  Input("input_complete.liquid");
  for(int istep=1; istep <= nstep; ++istep){
    Move();
    if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
    if(istep%10 == 0){
      Measure("liquid_out");
    }
  }
  Restart("liquid_out");
  FinalRestart("liquid_out");
  
  //SIMULAZIONE FASE GASSOSA
  
  Input("input_complete.gas");
  for(int istep=1; istep <= nstep; ++istep){
    Move();
    if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
    if(istep%10 == 0){
      Measure("gas_out");
    }
  }
  Restart("gas_out");
  FinalRestart("gas_out");
  
  return 0;
}

/***************************************************************/

void Input(std::string fileinput){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open(fileinput); //Read input

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
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;
  ReadInput >> nSimulations;
  ReadInput >> nStepfinal;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

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

/***************************************************************/

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

/***************************************************************/

void Move_WO_Vel(void){ //Move particles with Verlet algorithm
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
    
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];
    
    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

/***************************************************************/

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

/***************************************************************/

void Measure(std::string fileoutput){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open(fileoutput + "/" + "output_epot.dat",ios::app);
  Ekin.open(fileoutput + "/" + "output_ekin.dat",ios::app);
  Temp.open(fileoutput + "/" + "output_temp.dat",ios::app);
  Etot.open(fileoutput + "/" + "output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] );
     dy = Pbc( yold[i] - yold[j] );
     dz = Pbc( zold[i] - zold[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

/***************************************************************/

void MeasureR(std::string fileoutput, int i){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
  
  Epot.open(fileoutput + "/" + "output_epot" + to_string(i) + ".dat",ios::app);
  Ekin.open(fileoutput + "/" + "output_ekin" + to_string(i) + ".dat",ios::app);
  Temp.open(fileoutput + "/" + "output_temp" + to_string(i) + ".dat",ios::app);
  Etot.open(fileoutput + "/" + "output_etot" + to_string(i) + ".dat",ios::app);
  
  v = 0.0; //reset observables
  t = 0.0;
  
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
      
      dx = Pbc( xold[i] - xold[j] );
      dy = Pbc( yold[i] - yold[j] );
      dz = Pbc( zold[i] - zold[j] );
      
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      
      if(dr < rcut){
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        
        //Potential energy
        v += vij;
      }
    }
  }
  
  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  
  return;
}

/***************************************************************/

//Error function: statistical uncertainty estimation
double Error(std::vector<double> AV, std::vector<double> AV2, int n){
  if(n==0){
    return 0;
  }else{
    return sqrt((AV2[n]-pow(AV[n],2))/double(n));
  }
}

/***************************************************************/

//valuta incertezza con metodo a blocchi e stampa valori
void Eval_Print(std::string namefile, std::vector<double> AV, std::vector<double> AV2, int N){
  
  vector<int> x(N);
  iota(x.begin(), x.end(), 1);
  vector<double> sum_prog(N, 0.0);
  vector<double> su2_prog(N, 0.0);
  vector<double> err_prog(N, 0.0);
  
  for (int i=0; i<N; i++){
    for(int j=0; j<i+1; j++){
      sum_prog[i]+= AV[j];
      su2_prog[i]+= AV2[j];
    }
    sum_prog[i]/=(i+1);
    su2_prog[i]/=(i+1);
    err_prog[i]=Error(sum_prog, su2_prog, i );
  }
  
  ofstream Output;
  Output.open(namefile);
  if (Output.is_open()){
    for (int i=0; i<N; i++){
      Output << x[i] << "\t" << sum_prog[i]  << "\t" << err_prog[i] << endl;
    }
  } else cerr << "PROBLEM: Unable to open -> " << namefile << endl;
  Output.close();
}

/***************************************************************/

//valuta incertezza con metodo a blocchi e stampa valori in SI
void Eval_Print_SI(std::string namefile, std::vector<double> AV, std::vector<double> AV2, int N, bool energy){
  
  vector<int> x(N);
  iota(x.begin(), x.end(), 1);
  vector<double> sum_prog(N, 0.0);
  vector<double> su2_prog(N, 0.0);
  vector<double> err_prog(N, 0.0);
  
  for (int i=0; i<N; i++){
    for(int j=0; j<i+1; j++){
      sum_prog[i]+= AV[j];
      su2_prog[i]+= AV2[j];
    }
    sum_prog[i]/=(i+1);
    su2_prog[i]/=(i+1);
    err_prog[i]=Error(sum_prog, su2_prog, i );
  }
  
  ofstream Output;
  Output.open(namefile);
  if (Output.is_open()){
    if(energy==true){
      for (int i=0; i<N; i++){
        Output << x[i] << "\t" << sum_prog[i]*epsilon  << "\t" << err_prog[i]*epsilon << endl;
        }
      }else{
      for (int i=0; i<N; i++){
        Output << x[i] << "\t" << sum_prog[i]*epsilon/boltzmann  << "\t" << err_prog[i]*epsilon/boltzmann << endl;
      }
    }
  }else cerr << "PROBLEM: Unable to open -> " << namefile << endl;
  Output.close();
}

/***************************************************************/

void MeasureRB(std::string fileoutput, int i, int N, std::vector<std::vector<double> > &M)
{ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;
  
  //file valori istantanei
  Epot.open(fileoutput + "/" + "output_epot" + to_string(i) + ".dat",ios::app);
  Ekin.open(fileoutput + "/" + "output_ekin" + to_string(i) + ".dat",ios::app);
  Temp.open(fileoutput + "/" + "output_temp" + to_string(i) + ".dat",ios::app);
  Etot.open(fileoutput + "/" + "output_etot" + to_string(i) + ".dat",ios::app);
  //file valori medi
  
  v = 0.0; //reset observables
  t = 0.0;
  
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){
      
      dx = Pbc( xold[i] - xold[j] );
      dy = Pbc( yold[i] - yold[j] );
      dz = Pbc( zold[i] - zold[j] );
      
      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);
      
      if(dr < rcut){
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        
        //Potential energy
        v += vij;
      }
    }
  }
  
  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  
  M[0][N]+= stima_pot;
  M[2][N]+= stima_kin;
  M[4][N]+= stima_temp;
  M[6][N]+= stima_etot;
  
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  
  return;
}

/***************************************************************/

void ConfFinal(std::string fileoutput){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open(fileoutput + "/" +"config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

/***************************************************************/

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

/***************************************************************/

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

/***************************************************************/

void Restart(std::string fileoutput){
  
  if(nSimulations<1) {
      cout << "Error: nSimulations <1" << endl;
      abort();
  }
  
  if (restart==false) ConfFinal(fileoutput);
    else{
      
      int count=0;
      
      while(count<nSimulations-1){
      //Scrivo i risultati della simulazione precedente
        count++;
        cout << endl << endl << "Restart N. " << count << " of " << nSimulations << endl << endl;
        
        ofstream WriteConf;
        WriteConf.open(fileoutput + "/" + "old" + to_string(count) +".final");
        for (int i=0; i<npart; ++i){
          WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
        }
        WriteConf.close();
     
        WriteConf.open(fileoutput + "/" + "old" + to_string(count) +".0");
        for (int i=0; i<npart; ++i){
          WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
        }
        WriteConf.close();
        
        //calcolo passo r(t+dt);
        Move_WO_Vel();
        
        //calcolo v(t+dt/2)
        for(int i=0; i<npart; ++i){
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
        for(int i=0; i<npart; ++i){
          vx[i] = sqrt(c)*vx[i];
          vy[i] = sqrt(c)*vy[i];
          vz[i] = sqrt(c)*vz[i];
        }
        
        //stimo la nuova configurazione r(t)
        
        for(int i=0; i<npart; ++i){
          xold[i] = Pbc( x[i] - vx[i] * delta );
          yold[i] = Pbc( y[i] - vy[i] * delta );
          zold[i] = Pbc( z[i] - vz[i] * delta );
        }
        
        //avvio la simulazione con le nuove variabili
        int nconf = 1;
        for(int istep=1; istep <= nstep; ++istep){
          Move();           //Move particles with Verlet algorithm
          if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
          if(istep%10 == 0){
            MeasureR(fileoutput, count);     //Properties measurement
            //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
            nconf += 1;
          }
        }
      }
    }
  }

/***************************************************************/

void FinalRestart(std::string fileoutput){
  
  if (restart==false) ConfFinal(fileoutput);
  else{
    //Scrivo i risultati della simulazione precedente
  
    cout << endl << endl << "Restart N. " << nSimulations << " of " << nSimulations << endl << endl;
    
    ofstream WriteConf;
    WriteConf.open(fileoutput + "/" + "old" + to_string(nSimulations) +".final");
    for (int i=0; i<npart; ++i){
      WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    }
    WriteConf.close();
    
    WriteConf.open(fileoutput + "/" + "old" + to_string(nSimulations) +".0");
    for (int i=0; i<npart; ++i){
      WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    }
    WriteConf.close();
    
    //calcolo passo r(t+dt);
    Move_WO_Vel();
    
    //calcolo v(t+dt/2)
    for(int i=0; i<npart; ++i){
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
    for(int i=0; i<npart; ++i){
      vx[i] = sqrt(c)*vx[i];
      vy[i] = sqrt(c)*vy[i];
      vz[i] = sqrt(c)*vz[i];
    }
    
    //stimo la nuova configurazione r(t)
    
    for(int i=0; i<npart; ++i){
      xold[i] = Pbc( x[i] - vx[i] * delta );
      yold[i] = Pbc( y[i] - vy[i] * delta );
      zold[i] = Pbc( z[i] - vz[i] * delta );
    }
    
    //avvio la simulazione con le nuove variabili
    int     N=100, nconf=0;
    double  L=nStepfinal/(N);
    vector <vector<double> > M(8, vector<double>(N, 0.0));
    
    /*AVEpot
    AV2Epot
    AVEkin
    AV2Ekin
    AVTemp
    AV2Temp
    AVEtot
    AV2Etot*/
    
    for(int i=0; i<N ; i++){
      cout << "Number of time-steps: " << i*L << endl;
      for(int j=0; j<L; j++){
        Move();           //Move particles with Verlet algorithm
        if(j%10 == 0){ //misura i valori ogni 10 passi
          MeasureRB(fileoutput, nSimulations, i, M); //Properties measurement
          //        ConfXYZ(i+j/10);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
        }
      }
      //calcolo dei valori medi
      M[0][i]/=(L/10);
      M[1][i]=pow(M[0][i], 2);
      M[2][i]/=(L/10);
      M[3][i]=pow(M[2][i], 2);
      M[4][i]/=(L/10);
      M[5][i]=pow(M[4][i], 2);
      M[6][i]/=(L/10);
      M[7][i]=pow(M[6][i], 2);
    }
    Eval_Print_SI(fileoutput + "/" + "ave_epot.out", M[0], M[1], N, true);
    Eval_Print_SI(fileoutput + "/" + "ave_ekin.out", M[2], M[3], N, true);
    Eval_Print_SI(fileoutput + "/" + "ave_temp.out", M[4], M[5], N, false);
    Eval_Print_SI(fileoutput + "/" + "ave_etot.out", M[6], M[7], N, true);
    
    ConfFinal(fileoutput);
  }
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
