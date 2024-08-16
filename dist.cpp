#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include "dist.h"

using namespace std;

void dist(const int& total_n_atoms, vector<vector<double> >& box_dim, vector<vector<double> >& positions, vector<tuple<int, int, double, vector<PairwiseDistance> > >& pairwise_distances){
    //pairwise_distances.clear();
    for(int i=0; i<total_n_atoms; i++){
        for(int j=0; j<total_n_atoms; j++){
            //distances between the particles
            if(i!=j){
                double dx = positions[j][0] - positions[i][0];
                double dy = positions[j][1] - positions[i][1];
                double dz = positions[j][2] - positions[i][2];

                //PBC: minimum image convention
                dx -= box_dim[0][0] * round(dx / box_dim[0][0]);
                dy -= box_dim[1][1] * round(dy / box_dim[1][1]);
                dz -= box_dim[2][2] * round(dz / box_dim[2][2]);

                //distance r
                //double r = 0;
                double r = sqrt(dx*dx + dy*dy + dz*dz);
                if(r!=0){
                    vector<double> unit_r_vec = {dx/r, dy/r, dz/r};
                    PairwiseDistance pd = { i, j, r, unit_r_vec };
                    //all the pairs are listed out that have r < rc, repetition is there
                    pairwise_distances.push_back(make_tuple(i, j, r, vector<PairwiseDistance>{pd}));
        
                }
                
            }
        }
    }
}
