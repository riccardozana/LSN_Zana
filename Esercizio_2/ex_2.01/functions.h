#ifndef __Functions_h_
#define __Functions_h_

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>

using namespace std;

double Error(std::vector<double> , std::vector<double> , int );
void Eval_Print(std::string namefile, std::vector<double> AV, std::vector<double> AV2, int N, int L);

#endif
