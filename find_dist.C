#include <string>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


//positions from POSCAR
void find_dist(int& total_n_atoms, vector<vector<double>>& box_dim, vector<vector<double>>& positions, int& ith_particle, vector<double>& distances){
    //total number of atoms
    //L of the box_dim
    //positions
    //which particle is chosen
    
    //output the distance of ith particle wrt other particles in a 1D vector

    //lengths 

    for(int j=0; j<total_n_atoms; j++){
       
        double dx = 0, dy = 0, dz = 0;
        double r = 0;

        //distances between the particles
        dx = positions[j][0] - positions[ith_particle][0];
        dy = positions[j][1] - positions[ith_particle][1];
        dz = positions[j][2] - positions[ith_particle][2];

        //PBC: minimum image convention
        dx = dx - box_dim[0][0]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][0]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][0]*ceil(dz/box_dim[2][2]-0.5);
	    dy = dy - box_dim[0][1]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][1]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][1]*ceil(dz/box_dim[2][2]-0.5);
	    dz = dz - box_dim[0][2]*ceil(dx/box_dim[0][0]-0.5) - box_dim[1][2]*ceil(dy/box_dim[1][1]-0.5) - box_dim[2][2]*ceil(dz/box_dim[2][2]-0.5);
	
        //distance r
        r = sqrt(dx*dx + dy*dy + dz*dz);

        distances.push_back(r);
    }

}

