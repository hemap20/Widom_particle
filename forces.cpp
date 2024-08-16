//#include "pot_energy.h"
#include "forces.h"
#include "dist.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <tuple>
#include <map>

using namespace std;


double forces(const vector<tuple<int, int, double, vector<PairwiseDistance> > >& pairwise_distances, vector<tuple<int, double > >& pairwise_forces){
    //this method does not include repetition of pairs of particles. The value of the force calculated will be lesser than the actual value
    // Create a map to group j values by their corresponding i values
    map<int, vector<PairwiseDistance>> i_to_js_map;
    double F = 0.0;

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
        const vector<PairwiseDistance>& distances = pair.second;

        for (const auto& pd : distances) { //for every j that has the same i: within the rc radius
            F += 48 * (1.0 / pow(pd.r, 12) - 0.5 / pow(pd.r, 6));
        }   
    }
    return F;
}
