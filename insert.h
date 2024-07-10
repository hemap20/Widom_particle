#ifndef INSERT_H
#define INSERT_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tuple>
#include <gsl/gsl_rng.h>
#include "pairwise_dist.h"
#include "pot_energy.h"



using namespace std;


void insert_atom(gsl_rng * r,int& total_n_atoms, vector<vector<double>>& box_dim, vector<vector<double>>& positions);

#endif
