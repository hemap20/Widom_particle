#include "pe_i.h"
#include "pe_total.h"
#include <vector>
#include <iostream>

using namespace std;

double total_e (vector<vector<double> >& box_dim, const vector<vector<double> >& positions, const int N, double * vir) {
    int i;
    double pe_t = 0.0;

    for (i=0; i<N-1; i++) {
        for(int j=i+1; j<N; j++){
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
            //double s_12 = s * s * s * s * s * s * s * s * s * s * s * s;
            //double s_6 = s * s * s * s * s * s;
            pe_t += 4* (r6i * r6i - r6i);
            *vir += 48*(r6i*r6i-0.5*r6i);
        }
    }
    return pe_t;
}
