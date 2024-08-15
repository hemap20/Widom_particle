#ifndef MC_MOVE_H
#define MC_MOVE_H

#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include "pe_i.h"
#include "random.h"
#include "pe_total.h"
#include "mc_eq.h"

using namespace std;

void mc_move( double beta,vector<vector<double> > box_dim, vector<vector<double> >& positions,int total_n_atoms, double& step_size, int& trials, int& accepted_moves, int& total_accepted_moves, double& E, double& w, mt19937& gen);

#endif