#ifndef PAIRWISE_DIST_H
#define PAIRWISE_DIST_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

void dist(vector<vector<double>>& box_dim, vector<vector<double>>& positions,vector<tuple<int, int, double>>& pairwise_distances, int i, int j, double rc);

#endif
