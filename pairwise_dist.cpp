#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include "pairwise_dist.h"

using namespace std;


//compare with cutoff radius, check PBC, compute the distance
//start from i0
void dist(const int& total_n_atoms, double rc, vector<vector<double>>& box_dim, vector<vector<double>>& positions,vector<tuple<int, int, double, vector<PairwiseDistance>>>& pairwise_distances){
    cout << "pairwise_dist clear pw_dist.cpp" << endl; 
    pairwise_distances.clear();
    for(int i=0; i<total_n_atoms; i++){
        for(int j=0; j<total_n_atoms; j++){
            //distances between the particles
            if(i!=j){
                double dx = positions[j][0] - positions[i][0];
                double dy = positions[j][1] - positions[i][1];
                double dz = positions[j][2] - positions[i][2];

                //PBC: minimum image convention
                dx = dx - box_dim[0][0]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][0]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][0]*ceil(dz/box_dim[2][2]-0.5);
                dy = dy - box_dim[0][1]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][1]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][1]*ceil(dz/box_dim[2][2]-0.5);
                dz = dz - box_dim[0][2]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][2]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][2]*ceil(dz/box_dim[2][2]-0.5);
                
                //distance r
                double r = 0;
                r = sqrt(dx*dx + dy*dy + dz*dz);
                //cout << "r = " << r << endl;

                //cutoff radius
                if (r < rc) {
                    //cout << "pushing back pw_dist.cpp" << endl;
                    vector<double> unit_r_vec = {dx/r, dy/r, dz/r};
                    PairwiseDistance pd = { i, j, r, unit_r_vec };
                    //cout << " r vec: " << dx << " " << dy << " " << dz << endl; 
                    //all the pairs are listed out that have r < rc, repetition is there
                    pairwise_distances.push_back(make_tuple(i, j, r, vector<PairwiseDistance>{pd}));
                }
            }
        }
    }
}
