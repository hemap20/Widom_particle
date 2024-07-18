#include "pot_energy.h"
#include "pairwise_dist.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

// Only Potential Energy
//take values from pairwise_distances
double pot_energy(const vector<tuple<int, int, double, vector<PairwiseDistance>>>& pairwise_distances, double rc){
    
    double pot_energy_total = 0; 
    double rc_3 = 1.0 / (rc * rc * rc);
    double ecut = 4 * (rc_3 * rc_3 * rc_3 * rc_3 - rc_3 * rc_3);

    for (const auto& item : pairwise_distances) {
        double r = get<2>(item);
        double r2 = r * r;
        double r6i = 1.0 / (r2 * r2 * r2);
        pot_energy_total += 4 * (r6i * r6i - r6i) - ecut;
    }

    return pot_energy_total/2;
}