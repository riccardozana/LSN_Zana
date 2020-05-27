#ifndef __Functions_h_
#define __Functions_h_

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include "random.h"

using namespace std;

double Error(std::vector<double> , std::vector<double> , int );

void Eval_Print(std::string namefile, std::vector<double> AV, std::vector<double> AV2, int N, int L);

double S_direct(double S_0, double mu, double sgm, double T, Random * rnd);

double S_discret(double S_0, double mu, double sgm, double T, Random * rnd, int N);

double C00(double S_T, double r, double T,  double K);

double P00(double S_T, double r, double T,  double K);

#endif
