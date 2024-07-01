#include "pot_energy.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

struct PairwiseDistance {
    int i;
    int j;
    double r;
    vector<double> unit_r_vec;
};

// Only Potential Energy
double pot_energy(int N, double rc, vector<vector<double>>& box_dim, vector<vector<double>>& positions,
                        vector<tuple<int, int, double, vector<PairwiseDistance>>>& pairwise_distances, double rho) {
   
    double pot_energy_total = 0.0;
    double rc_3 = 1.0 / (rc * rc * rc);
    double ecut = 4 * (rc_3 * rc_3 * rc_3 * rc_3 - rc_3 * rc_3);

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            // Distance between particle i and j
            double dx = positions[j][0] - positions[i][0];
            double dy = positions[j][1] - positions[i][1];
            double dz = positions[j][2] - positions[i][2];

            // Periodic boundary conditions
            dx = dx - box_dim[0][0]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][0]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][0]*ceil(dz/box_dim[2][2]-0.5);
            dy = dy - box_dim[0][1]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][1]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][1]*ceil(dz/box_dim[2][2]-0.5);
            dz = dz - box_dim[0][2]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][2]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][2]*ceil(dz/box_dim[2][2]-0.5);

            // Distance r
            double r = sqrt(dx * dx + dy * dy + dz * dz);

            // Unit vector
            vector<double> unit_r_vec = { dx / r, dy / r, dz / r };

            // Within cutoff radius
            if (r < rc) {
                PairwiseDistance pd = { i, j, r, unit_r_vec };
                pairwise_distances.push_back(make_tuple(i, j, r, vector<PairwiseDistance>{pd}));

                double r2 = r * r;
                double r6i = 1.0 / (r2 * r2 * r2);
                pot_energy_total += 4 * (r6i * r6i - r6i) - ecut;
            }
        }
    }

    return pot_energy_total;
}

