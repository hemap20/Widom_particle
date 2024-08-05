#ifndef RAND_POS_H
#define RAND_POS_H

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

Coordinates pos(vector<vector<double> >& box_dim, double step_size);

#endif