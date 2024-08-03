#ifndef POSITIONS_H
#define POSITIONS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "positions.h"


using namespace std;

struct Particle {
    double x, y, z;
};

vector<Particle> generateParticles(int n, double density);

#endif