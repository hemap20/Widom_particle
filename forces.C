#include <cmath>
#include <iostream>
#include <vector>
#include <tuple>
#include "pot_energy.h"

using namespace std;

struct PairwiseForce {
    int i;
    //int j;
    //double F;
    vector<double> F_vec;
};


void force(const vector<PairwiseDistance>& pairwise_distances,vector<tuple<int, int, double, vector<PairwiseForce>>>& pairwise_forces, int N ){
    double F = 0;
    // Calculate forces based on pairwise distances
    pairwise_forces.clear();
    for (const auto& tuple : pairwise_distances) {
        int i = tuple.i;
        int j = tuple.j;
        double r = tuple.r;
        const vector<double>& unit_r_vec = tuple.unit_r_vec;

        double Fx_i = 0, Fy_i = 0, Fz_i = 0;
        //double Fx_j = 0, Fy_j = 0, Fz_j = 0;

        for (int k = 0; k<pairwise_distances.size(); k++) {
            // Calculate force magnitude
            double F = 48 * (1.0 / pow(pd.r, 8) - 0.5 / pow(pd.r, 4));

            // Apply force components alongdouble force(const vector<PairwiseDistance>& pairwise_distances,vector<tuple<int, int, double, vector<PairwiseForce>>>& pairwise_forces, int N ); each axis
            // Force on particle i
            Fx_i += pd.unit_r_vec[0] * F; // Force component along x-direction
            Fy_i += pd.unit_r_vec[1] * F; // Force component along y-direction
            Fz_i += pd.unit_r_vec[2] * F; // Force component along z-direction

            // Accumulate forces on particle i
            // Fx[i] += Fx_i;
            // Fy[i] += Fy_i;
            // Fz[i] += Fz_i;

            // Force on particle j (opposite direction)
            // Fx_j += -pd.unit_r_vec[0] * F; // Force component along x-direction for particle j
            // Fy_j += -pd.unit_r_vec[1] * F; // Force component along y-direction for particle j
            // Fz_j += -pd.unit_r_vec[2] * F; // Force component along z-direction for particle j

            // Accumulate forces on particle j 
            // Fx[j] += Fx_j;
            // Fy[j] += Fy_j;
            // Fz[j] += Fz_j;
        }
        vector<double> F_vector = { Fx_i, Fy_i, Fz_i};
        //PairwiseForce pd = { i, F_vector }; //stores the value of index i, total force vector that i experiences 
        //PairwiseForce pd = { i, F_vector };
        pairwise_forces.push_back(make_tuple(i, F_vector));

        
    }
}





