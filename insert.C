#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tuple>
#include <gsl/gsl_rng.h>
//#include </home/hema/Codes/basics/Basics/gsl-2.8/gsl/gsl_rng.h>
#include "pairwise_dist.h"
#include "pot_energy.h"
#include "insert.h"



using namespace std;

//potential energy for an entire system of particles
      //  generate a random set of coordinates from r
      //  add to the positions list, update total_number_atoms
void insert_atom(gsl_rng * r,int& total_n_atoms, vector<vector<double>>& box_dim, vector<vector<double>>& positions){
    double rx = (gsl_rng_uniform(r)-0.5)*box_dim[0][0];
    double ry = (gsl_rng_uniform(r)-0.5)*box_dim[1][1];
    double rz = (gsl_rng_uniform(r)-0.5)*box_dim[2][2];

    positions.push_back({ rx, ry, rz });
    total_n_atoms = positions.size();

}
