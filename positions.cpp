#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>
#include "positions.h"

using namespace std;


void generateParticles(vector<vector<double> > & positions, int n, double density, vector<vector<double> >& box_dim) {
    // Calculate the volume of the system
    double volume = n / density;
    cout << "vol" << volume << endl;

    // Assuming a cubic volume for simplicity
    double Lx = box_dim[0][0];
    double Ly = box_dim[1][1];
    double Lz = box_dim[2][2];

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis_x(0.0, Lx);
    uniform_real_distribution<> dis_y(0.0, Ly);
    uniform_real_distribution<> dis_z(0.0, Lz);

    for (int i = 0; i < n; ++i) {
        positions[i][0] = dis_x(gen);
        positions[i][1] = dis_y(gen);
        positions[i][2] = dis_z(gen);
    }

}

