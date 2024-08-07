#ifndef MC_EQ_H
#define MC_EQ_H

#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include "pe_i.h"
#include "random.h"

using namespace std;

void mc_eq( double e, double beta, double s, vector<vector<double> > box_dim, vector<vector<double> >& positions,int total_n_atoms,double& step_size, int& trials, int& accepted_moves, int& total_accepted_moves, int seed, double& E, double& sum_E);

#endif