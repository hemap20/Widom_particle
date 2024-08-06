#ifndef MC_EQ_H
#define MC_EQ_H

#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include "pe_i.h"
#include "random.h"

using namespace std;

void mc_eq(const double e,const double beta,const double s,const vector<vector<double> > box_dim, vector<vector<double> >& positions,int total_n_atoms, double& w, double& step_size, int& trials, int& accepted_moves,  int seed, double& E);

#endif