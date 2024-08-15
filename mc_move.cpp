#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include "pe_i.h"
#include "pe_total.h"
#include "random.h"
#include "mc_eq.h"
#include "mc_move.h"

using namespace std;

void mc_move( double beta, vector<vector<double> > box_dim, vector<vector<double> >& positions,int total_n_atoms, double& step_size, int& trials, int& accepted_moves, int& total_accepted_moves, double& E, double& w, mt19937& gen){
    int i = rand() % total_n_atoms;

    // Calculate the initial energy of atom i
    double en_0 = e_i(0, box_dim, positions, i, total_n_atoms);

    // Generate a newpositions position for atom i
    Coordinates new_pos = random(box_dim, step_size,gen);
    double x_new = new_pos.x;
    double y_new = new_pos.y;
    double z_new = new_pos.z;        
    
    positions.push_back({x_new, y_new, z_new});
    total_n_atoms = total_n_atoms + 1;

    // Calculate the new energy of atom i
    double en_new = e_i(0,box_dim,positions, total_n_atoms-1, total_n_atoms);

    uniform_real_distribution<> dis_real(0.0, 1.0);
    //mt19937 gen(seed);
    double R = dis_real(gen);

    // Metropolis acceptance criterion
    if ((en_new < en_0)||(R < exp(-beta * (en_new - en_0)))) {
        w += exp(-beta * en_new);
        //cout << "w " << w << endl;
        accepted_moves++; //register the insertion
        total_accepted_moves++;
        E += en_new - en_0;
    }
    positions.pop_back();
    total_n_atoms = total_n_atoms - 1;

    trials++;
    
    double acceptance_ratio = static_cast<double>(accepted_moves) / trials;
    
    if (acceptance_ratio < 0.3) {
        step_size *= 0.3; // Decrease step size
        // Reset counters
        trials = 0;
        accepted_moves = 0;
        //cout << "Step: " << total_trials << ", Acceptance Ratio: " << acceptance_ratio << endl;
    } else if (acceptance_ratio > 0.5) {
        step_size *= 1.2; // Increase step size
        // Reset counters
        trials = 0;
        accepted_moves = 0;
        //cout << "Step: " << total_trials << ", Acceptance Ratio: " << acceptance_ratio << endl;
    }    

}