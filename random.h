#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <set>

using namespace std;

struct Coordinates {
    double x, y, z;
};

extern set<double> generatedDoubles;

Coordinates random(vector<vector<double> >& box_dim, double step_size, mt19937& gen);

#endif