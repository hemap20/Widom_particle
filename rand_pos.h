#ifndef RAND_POS_H
#define RAND_POS_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tuple>
#include <random>
#include <set>

using namespace std;

set<double> generatedDoubles;

tuple<double, double, double> pos(int& total_n_atoms, vector<vector<double> >& box_dim, vector<vector<double> >& positions, double step_size);

#endif 