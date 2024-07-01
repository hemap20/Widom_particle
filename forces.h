#ifndef FORCES_H
#define FORCES_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
//include path to pairwise_distance


using namespace std; 

struct PairwiseForce {
    int i;
    int j;
    double F;
    vector<double> F_vec;
};

double force(const vector<PairwiseDistance>& pairwise_distances,vector<tuple<int, int, double, vector<PairwiseForce>>>& pairwise_forces, int N );

#endif 