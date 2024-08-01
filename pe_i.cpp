#include "pe_i.h"
#include <iostream>
#include <vector>

using namespace std;

double e_i (int i0, const double e, vector<vector<double>>& box_dim, const double s, const vector<vector<double>>& positions, int i, const int N) {
    int j = 0;
    double pe_i = 0;
    
    for(j=i0; j<N; j++){
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
            double r2 = r * r;
            double r6i = 1.0 / (r2 * r2 * r2);
            double s_12 = s * s * s * s * s * s * s * s * s * s * s * s;
            double s_6 = s * s * s * s * s * s;
            pe_i += 4*e* (s_12*r6i * r6i - s_6*r6i);
        }
    }
    return pe_i;
}