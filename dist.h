#ifndef DIST_H
#define DIST_H

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

// Corrected declaration of the dist function
void dist(const int& total_n_atoms, vector<vector<double>>& box_dim, vector<vector<double>>& positions,vector<tuple<int, int, double, vector<PairwiseDistance>>>& pairwise_distances);

#endif
