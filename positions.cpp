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

vector<Particle> generateParticles(int n, double density) {
    // Calculate the volume of the system
    double volume = n / density;

    // Assuming a cubic volume for simplicity
    double side_length = cbrt(volume);

    vector<Particle> particles;
    particles.reserve(n);

    // Seed the random number generator
    srand(time(nullptr));

    for (int i = 0; i < n; ++i) {
        Particle p;
        p.x = static_cast<double>(std::rand()) / RAND_MAX * side_length;
        p.y = static_cast<double>(std::rand()) / RAND_MAX * side_length;
        p.z = static_cast<double>(std::rand()) / RAND_MAX * side_length;
        particles.push_back(p);
    }

    return particles;
}
