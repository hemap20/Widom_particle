#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>
#include <set>

using namespace std;

void trialmove(const double beta, const double e, const double s, vector<vector<double> >& box_dim, vector<vector<double>>& positions, int total_n_atoms) {
    // Select a random atom
    int i = rand() % total_n_atoms;

    // Calculate the initial energy of atom i
    double en_0 = e_i(0, e, box_dim, s, positions, i, total_n_atoms);

    double x_old = positions[i][0];
    double y_old = positions[i][1];
    double z_old = positions[i][2];

    double x_new = 0, y_new = 0, z_new = 0;
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
            x_new = abs(coord);
            if (x_new < 0.001) x_new = 0.001;
            if (x_new > box_dim[0][0]) x_new = box_dim[0][0];
        }
        else if(i==1){
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_y(0.001, box_dim[1][1]);
                newDouble = dis_y(gen);
                coord = newDouble + step_size * (dis_y(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            y_new = abs(coord);
            if (y_new < 0.001) ry = 0.001;
            if (y_new > box_dim[1][1]) y_new = box_dim[1][1];
        }
        else{
            double newDouble;
            double coord;
            do {
                static uniform_real_distribution<> dis_z(0.001, box_dim[2][2]);
                newDouble = dis_z(gen);
                coord = newDouble + step_size * (dis_z(gen) - 0.5);
            } while (!generatedDoubles.insert(newDouble).second);
            z_new = abs(coord);
            if (z_new < 0.001) rz = 0.001;
            if (z_new > box_dim[2][2]) z_new = box_dim[2][2];
        }
    }
    positions[i][0] = x_new;
    positions[i][1] = y_new;
    positions[i][2] = z_new;

    // Calculate the new energy of atom i
    double en_new = e_i(0, e, box_dim, s, positions, i, total_n_atoms);

    uniform_real_distribution<> dis_real(0.0, 1.0);
    double R = dis_real(gen);

    // Metropolis acceptance criterion
    if (R < exp(-beta * (en_new - en_0))) {
        // Accept the move, position is already updated

    } else {
        // Reject the move, revert to the original position
        positions[i][0] = x_old;
        positions[i][1] = y_old;
        positions[i][2] = z_old;
    }
}
