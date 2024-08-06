#include <iostream>
#include <vector>
#include <random>
#include <set>
#include "random.h"

using namespace std;

set<double> generatedDoubles;

Coordinates random(vector<vector<double> >& box_dim, double step_size){

    double x = 0, y = 0, z = 0;
    static random_device rd;  // Initialize random device
    static mt19937 gen(rd()); // Seed the generator
    
    for(int j=0; j<3; j++){
        if(j==0){
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_x(0.001, box_dim[0][0]);
                newDouble = dis_x(gen);
                coord = newDouble + step_size * (dis_x(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            x = abs(coord);
            if (x < 0.001) x = 0.001;
            if (x > box_dim[0][0]) x = box_dim[0][0];
        }
        else if(j==1){
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_y(0.001, box_dim[1][1]);
                newDouble = dis_y(gen);
                coord = newDouble + step_size * (dis_y(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            y = abs(coord);
            if (y < 0.001) y = 0.001;
            if (y > box_dim[1][1]) y = box_dim[1][1];
        }
        else{
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_z(0.001, box_dim[2][2]);
                newDouble = dis_z(gen);
                coord = newDouble + step_size * (dis_z(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            z = abs(coord);
            if (z < 0.001) z = 0.001;
            if (z > box_dim[2][2]) z = box_dim[2][2];
        }
    }
    return {x, y, z};

}