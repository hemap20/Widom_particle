#include "pairwise_dist.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

struct PairwiseDistance {
    int i;
    int j;
    double r;
    vector<double> unit_r_vec;
};

//positions from POSCAR

//compare with cutoff radius, check PBC, compute the distance
void dist(int N, double rc, vector<vector<double>>& box_dim, vector<vector<double>>& positions,vector<tuple<int, int, double, vector<PairwiseDistance>>>& pairwise_distances, double rho){
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            //distances between the particles
            double dx = positions[j][0] - positions[i][0];
            double dy = positions[j][1] - positions[i][1];
            double dz = positions[j][2] - positions[i][2];

            //PBC: minimum image convention
            dx = dx - box_dim[0][0]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][0]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][0]*ceil(dz/box_dim[2][2]-0.5);
            dy = dy - box_dim[0][1]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][1]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][1]*ceil(dz/box_dim[2][2]-0.5);
            dz = dz - box_dim[0][2]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][2]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][2]*ceil(dz/box_dim[2][2]-0.5);

            //distance r
            double r = sqrt(dx*dx + dy*dy + dz*dz);

            vector<double> unit_r_vec = {dx/r, dy/r, dz/r};

            //cutoff radius
            if (r < rc) {
                PairwiseDistance pd = { i, j, r, unit_r_vec };
                pairwise_distances.push_back(make_tuple(i, j, r, vector<PairwiseDistance>{pd}));
            }

        }
    }

    
}

