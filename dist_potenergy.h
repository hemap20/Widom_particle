#ifndef DIST_POTENERGY_H
#define DIST_POTENERGY_H

#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

double dist_potenergy(int N, double rc, vector<double> &Fx, vector<double> &Fy, vector<double> &Fz,
                vector<vector<double>>& box_dim, vector<vector<double>>& positions,vector<tuple<int, int, double>>& pairwise_distances, double rho);

#endif
