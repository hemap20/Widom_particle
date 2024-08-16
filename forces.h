#ifndef FORCES_H
#define FORCES_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <map>
//include path to pairwise_distance


using namespace std; 

struct PairwiseDistance;

struct PairwiseForce {
    int i;
    double F;
    vector<double> F_vec;
}; 

void forces(const vector<tuple<int, int, double, vector<PairwiseDistance> > >& pairwise_distances, vector<tuple<int, double, vector<PairwiseForce> > >& pairwise_forces);

#endif 