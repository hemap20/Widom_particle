#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tuple>
#include <random>
#include <set>
//#include </home/hema/Codes/basics/Basics/gsl-2.8/gsl/gsl_rng.h>
#include "pairwise_dist.h"
#include "pot_energy.h"
#include "insert.h"

using namespace std;

set<double> generatedDoubles;

//potential energy for an entire system of particles
      //  generate a random set of coordinates from r
      //  add to the positions list, update total_number_atoms
void insert_atom(int& total_n_atoms, vector<vector<double> >& box_dim, vector<vector<double> >& positions){
    double rx = 0, ry = 0, rz = 0;
    
    static uniform_real_distribution<> dis(0.0, 1.0);
    static mt19937 gen(std::random_device{}());
    
    for(int i=0; i<3; i++){
        if(i==0){
            double newDouble;
            do {
                newDouble = dis(gen);
            } while (!generatedDoubles.insert(newDouble).second);
            rx = newDouble;
        }
        else if(i==1){
            double newDouble;
            do {
                newDouble = dis(gen);
            } while (!generatedDoubles.insert(newDouble).second);
            ry = newDouble;
        }
        else{
            double newDouble;
            do {
                newDouble = dis(gen);
            } while (!generatedDoubles.insert(newDouble).second);
            rz = newDouble;
        }
    }
   
    rx = abs((rx - 0.5)*box_dim[0][0]);
    ry = abs((ry - 0.5)*box_dim[1][1]);
    rz = abs((rz - 0.5)*box_dim[2][2]);
    
    cout << "pushing back insert.cpp" << endl;
    positions.push_back({ rx, ry, rz });
    total_n_atoms = positions.size();

}