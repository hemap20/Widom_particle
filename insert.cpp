#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tuple>
#include <random>
#include <set>
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
    
    static mt19937 gen(random_device{}());
    
    for(int i=0; i<3; i++){
        if(i==0){
            double newDouble;
            do {
                static uniform_real_distribution<> dis_x(0.001, box_dim[0][0]);
                newDouble = dis_x(gen);
            } while (!generatedDoubles.insert(newDouble).second);
            rx = abs(newDouble);
        }
        else if(i==1){
            double newDouble;
            do {
                static uniform_real_distribution<> dis_y(0.001, box_dim[1][1]);
                newDouble = dis_y(gen);
            } while (!generatedDoubles.insert(newDouble).second);
            ry = abs(newDouble);
        }
        else{
            double newDouble;
            do {
                static uniform_real_distribution<> dis_z(0.001, box_dim[2][2]);
                newDouble = dis_z(gen);
            } while (!generatedDoubles.insert(newDouble).second);
            rz = abs(newDouble);
        }
    }
   
    positions.push_back({ rx, ry, rz });
    total_n_atoms = positions.size();

}