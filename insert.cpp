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

void insert_atom(int& total_n_atoms, vector<vector<double> >& box_dim, vector<vector<double> >& positions, double step_size){
    double rx = 0, ry = 0, rz = 0;
    
    static mt19937 gen(random_device{}());
    
    for(int i=0; i<3; i++){
        if(i==0){
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_x(0.001, box_dim[0][0]);
                newDouble = dis_x(gen);
                coord = newDouble + step_size * (dis_x(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            rx = abs(coord);
            if (rx < 0.001) rx = 0.001;
            if (rx > box_dim[0][0]) rx = box_dim[0][0];
        }
        else if(i==1){
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_y(0.001, box_dim[1][1]);
                newDouble = dis_y(gen);
                coord = newDouble + step_size * (dis_y(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            ry = abs(coord);
            if (ry < 0.001) ry = 0.001;
            if (ry > box_dim[1][1]) ry = box_dim[1][1];
        }
        else{
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_z(0.001, box_dim[2][2]);
                newDouble = dis_z(gen);
                coord = newDouble + step_size * (dis_z(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            rz = abs(coord);
            if (rz < 0.001) rz = 0.001;
            if (rz > box_dim[2][2]) rz = box_dim[2][2];
        }
    }
   
    positions.push_back({ rx, ry, rz });
    total_n_atoms = positions.size();

}