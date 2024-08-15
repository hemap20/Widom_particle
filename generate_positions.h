#ifndef GENERATE_POSITIONS_H
#define GENERATE_POSITIONS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <random>
#include <ctime>

using namespace std;

void generate_positions(vector<vector<double> > & positions, int n, vector<vector<double> >& box_dim, mt19937& gen);

#endif