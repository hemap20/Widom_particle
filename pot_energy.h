#ifndef POT_ENERGY_H
#define POT_ENERGY_H

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

double pot_energy(int N, double rc,vector<vector<double>>& box_dim, vector<vector<double>>& positions,
                        vector<tuple<int, int, double, vector<PairwiseDistance>& pairwise_distances, double rho);

#endif 