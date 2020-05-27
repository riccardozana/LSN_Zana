#ifndef __Functions_h_
#define __Functions_h_

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include "random.h"
#include "funzioneBase.h"

using namespace std;

double Error(std::vector<double> , std::vector<double> , int );

void Eval_Print(std::string namefile, std::vector<double> AV, std::vector<double> AV2, int N);
void Eval_Store(std::vector<double> &v, std::vector<double> AV, std::vector<double> AV2, int N, int L);

double R2(std::vector<double> v);

void T_unif (vector<double> & xold, double radius, Random* rnd);
void T_gauss (vector<double> & xold, double sigma, Random* rnd);

double Eval_Alpha(vector<double> xold, vector<double> x, FunzioneBase2* p);

void Accept (double alpha, vector<double> & xold, vector<double> x, int & count, Random * rnd);

void Equilibrate_unif_print (vector<double> &xold, FunzioneBase2 * f, double radius, int N, Random * rnd);
void Equilibrate_unif (vector<double> &xold, FunzioneBase2 * f, double radius, int N, Random * rnd);
void Equilibrate_gauss (vector<double> &xold, FunzioneBase2 * f, double sigma, int N, Random * rnd);

#endif
