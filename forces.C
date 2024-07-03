#include "pot_energy.h"
#include "forces.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <tuple>
#include <map>

using namespace std;

// struct PairwiseForce {
//     int i;
//     vector<double> F_vec;
// };

void forces(const vector<tuple<int, int, double, vector<PairwiseDistance>>>& pairwise_distances, vector<tuple<int, double, vector<PairwiseForce>>>& pairwise_forces) {
    //this method does not include repetition of pairs of particles. The value of the force calculated will be lesser than the actual value
    // Create a map to group j values by their corresponding i values
    map<int, vector<PairwiseDistance>> i_to_js_map;

    // Populate the map
    for (const auto& tuple : pairwise_distances) {
        const vector<PairwiseDistance>& distances = get<3>(tuple);
        for (const auto& pd : distances) {
            i_to_js_map[pd.i].push_back(pd);
        }
    }

    // Iterate through the map to loop over unique i values and their corresponding j values
    for (const auto& pair : i_to_js_map) { //for every unique i
        int i = pair.first;
        double Fx_i = 0, Fy_i = 0, Fz_i = 0, F = 0;
        const vector<PairwiseDistance>& distances = pair.second;

        for (const auto& pd : distances) { //for every j that has the same i: within the rc radius
            F = 48 * (1.0 / pow(pd.r, 8) - 0.5 / pow(pd.r, 4));
            Fx_i += pd.unit_r_vec[0] * F; // Force component along x-direction
            Fy_i += pd.unit_r_vec[1] * F; // Force component along y-direction
            Fz_i += pd.unit_r_vec[2] * F; // Force component along z-direction
            //summing the forces for each i in the list
        }
        vector<double> F_vector = { Fx_i, Fy_i, Fz_i};
        PairwiseForce pf = { i, F, F_vector };
        pairwise_forces.push_back(make_tuple(i, F, vector<PairwiseForce>{pf}));
        //pairwise_forces.push_back(make_tuple(i, vector<PairwiseForce>{PairwiseForce{i, F_vector}}));
    }
}





