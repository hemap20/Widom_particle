#include "dist_potenergy.h"
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;


//positions from POSCAR

//compare with cutoff radius, check PBC, compute the distance
double dist_potenergy(int N, double rc, vector<double> &Fx, vector<double> &Fy, vector<double> &Fz,
                vector<vector<double>>& box_dim, vector<vector<double>>& positions,vector<tuple<int, int, double>>& pairwise_distances, double rho){
   
    double dx = 0, dy = 0, dz = 0;
    double r = 0;
    double pot_energy = 0.0, F = 0.0;
    double r2, r6i, rc_3, ecorr, ecut;

    //zeroing the forces
    for(int i=0; i<N; i++){
        Fx[i] = Fy[i] = Fz[i] = 0.0;
    }
   
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            //distance between particle i and j
            double dx = 0, dy = 0, dz = 0;
            double r = 0;

            dx = positions[j][0] - positions[i][0];
            dy = positions[j][1] - positions[i][1];
            dz = positions[j][2] - positions[i][2];

            //Periodic boundary conditions: applying minimum image convention
            dx = dx - box_dim[0][0]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][0]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][0]*ceil(dz/box_dim[2][2]-0.5);
            dy = dy - box_dim[0][1]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][1]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][1]*ceil(dz/box_dim[2][2]-0.5);
            dz = dz - box_dim[0][2]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][2]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][2]*ceil(dz/box_dim[2][2]-0.5);

            //distance r
            r = sqrt(dx*dx + dy*dy + dz*dz);
            
            //square of the distance between the particles
            r2 = r*r;

            //within the cutoff radius
            if(r<rc){
                pairwise_distances.push_back(make_tuple(i, j, r));
                rc_3 = 1.0/(rc*rc*rc);
                ecorr = 8*3.145*rho*(rc_3*rc_3*rc_3/9.0-rc_3/3.0);
                ecut = 4*(rc_3*rc_3*rc_3*rc_3-rc_3*rc_3);

                
                r6i         = 1.0/(r2*r2*r2);
                pot_energy  += 4*(r6i*r6i - r6i) - ecut; //e_cut is subtracted to keep the graph of energy continous at the cutoff radius distance   
                F           = 48*(r6i*r6i-0.5*r6i); //force between the i and j particles, not the total force

                //calculating the components of the forces and adding them up 
                Fx[i] += dx*F/r2;
                Fx[j] -= dx*F/r2;
                Fy[i] += dy*F/r2;
                Fy[j] -= dy*F/r2;
                Fz[i] += dz*F/r2;
                Fz[j] -= dz*F/r2;
            }
        }
    }
    return pot_energy+ N*ecorr; //potential energy 
}







