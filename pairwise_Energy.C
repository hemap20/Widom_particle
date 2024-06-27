#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <tuple>
#include <vector>

//coordinates: ri, forces: fi

//use formula to calculate ecorr and ecut

double potential_energy(double& Fx, double& Fy, double& Fz, vector<tuple<int, int, double>>& pairwise_distances, 
                        int N, double L, double rc2, double e_corr, double e_cut, double *vir){
    
    double pot_energy = 0.0, F = 0.0;
    double r2, r6i;

    //zeroing the forces
    for(int i=0; i<N; i++){
        Fx[i] = Fy[i] = Fz[i] = 0.0;
    }
    *vir = 0.0;

    for (const auto& distance_info : pairwise_distances) {
        //cout << "i: " << get<0>(distance_info) << ", j: " << get<1>(distance_info) << ", r: " << get<2>(distance_info) << endl;
        int i = get<0>(distance_info);
        int j = get<1>(distance_info);
        double r = get<2>(distance_info);
        r2 = r*r;
        r6i = 1.0/(r2*r2*r2);
        pot_energy += 4*(r6i*r6i - r6i) - e_cut; //e_cut is subtracted to keep the graph of energy continous at the cutoff radius distance   
        F = 48*(r6i*r6i-0.5*r6i);

        //calculating the components of the forces and adding them up
        Fx[i] += dx*F/r2;
        Fx[j] -= dx*F/r2;
        Fy[i] += dy*F/r2;
        Fy[j] -= dy*F/r2;
        Fz[i] += dz*F/r2;
        Fz[j] -= dz*F/r2;

        *vir += F; 

    }
    return pot_energy + N*e_cor; //ecor = energy due to tail corrections 
}

    // for(int i=0; i<N; i++){
    //     for(int j=i+1; j<N; j++){
    //         //distance between particle i and j
    //         dx = Rx[i] - Rx[j];
    //         dy = Ry[i] - Ry[j];
    //         dz = Rz[i] - Rz[j];

    //         //Periodic boundary conditions: applying minimum image convention
    //         if(dx>hL)       dx-=L;
    //         else if (dx<hL) dx+=L;
    //         if(dy>hL)       dy-=L;
    //         else if(dy<hL)  dy+=L;
    //         if(dz>hL)       dz-=L;
    //         else if(dz<hL)  dz+=L;  

    //         //square of the distance between the particles
    //         r2 = dx*dx + dy*dy + dz*dz;

    //         //within the cutoff radius
    //         if(r2<rc2){
    //             r6i         = 1.0/(r2*r2*r2);
    //             pot_energy  += 4*(r6i*r6i - r6i) - e_cut; //e_cut is subtracted to keep the graph of energy continous at the cutoff radius distance   
    //             F           = 48*(r6i*r6i-0.5*r6i);

    //             //calculating the components of the forces and adding them up
    //             Fx[i] += dx*F/r2;
    //             Fx[j] -= dx*F/r2;
    //             Fy[i] += dy*F/r2;
    //             Fy[j] -= dy*F/r2;
    //             Fz[i] += dz*F/r2;
    //             Fz[j] -= dz*F/r2;

    //             //adding up the internal particle forces
    //             *vir += F; 
    //         }
    //     }
    // }
    // return pot_energy + N*e_cut; //potential energy 


