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

//potential energy for an entire system of particles
      //  generate a random set of coordinates from r
      //  add to the positions list, update total_number_atoms
void insert_atom(int seed, int& total_n_atoms, vector<vector<double> >& box_dim, vector<vector<double> >& positions){
    double rx = 0, ry = 0, rz = 0;
    
    uniform_real_distribution<> dis_x(0.0, 1.0);
    mt19937 gen_x(seed);
    rx = dis_x(gen_x);

    uniform_real_distribution<> dis_y(0.0, 1.0);
    mt19937 gen_y(seed+1);
    ry = dis_y(gen_y);

    uniform_real_distribution<> dis_z(0.0, 1.0);
    mt19937 gen_z(seed+2);
    rz = dis_z(gen_z);

    rx = abs((rx - 0.5)*box_dim[0][0]);
    ry = abs((ry - 0.5)*box_dim[1][1]);
    rz = abs((rz - 0.5)*box_dim[2][2]);
    
    cout << "pushing back insert.cpp" << endl;
    positions.push_back({ rx, ry, rz });
    total_n_atoms = positions.size();

}