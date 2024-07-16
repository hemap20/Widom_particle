#ifndef INSERT_H
#define INSERT_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tuple>
#include <random>
//#include </home/hema/Codes/basics/Basics/gsl-2.8/gsl/gsl_rng.h>
#include "pairwise_dist.h"
#include "pot_energy.h"
#include "insert.h"

using namespace std;

void insert_atom(int seed, int& total_n_atoms, vector<vector<double> >& box_dim, vector<vector<double> >& positions);

#endif
