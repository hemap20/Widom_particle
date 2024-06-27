#include "pairwise_dist.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;


//positions from POSCAR

//compare with cutoff radius, check PBC, compute the distance
void dist(vector<vector<double>>& box_dim, vector<vector<double>>& positions,vector<tuple<int, int, double>>& pairwise_distances, int i, int j, double rc){
    double dx = 0, dy = 0, dz = 0;
    double r = 0;

    vector<double> vec1D;
    vector<double> vec2D;

    //distances between the particles
    dx = positions[j][0] - positions[i][0];
    dy = positions[j][1] - positions[i][1];
    dz = positions[j][2] - positions[i][2];

    //PBC: minimum image convention
    dx = dx - box_dim[0][0]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][0]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][0]*ceil(dz/box_dim[2][2]-0.5);
    dy = dy - box_dim[0][1]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][1]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][1]*ceil(dz/box_dim[2][2]-0.5);
    dz = dz - box_dim[0][2]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][2]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][2]*ceil(dz/box_dim[2][2]-0.5);

    //distance r
    r = sqrt(dx*dx + dy*dy + dz*dz);

    //cutoff radius
    if (r < rc) {
        pairwise_distances.push_back(make_tuple(i, j, r));
    }
}

