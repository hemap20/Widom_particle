#ifndef POT_ENERGY_H
#define POT_ENERGY_H

#include "pot_energy.h"
#include "pairwise_dist.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

double pot_energy(const double e, const double s, const vector<tuple<int, int, double, vector<PairwiseDistance> > >& pairwise_distances);

#endif 