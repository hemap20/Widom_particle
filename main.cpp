#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include <random>
#include "pe_i.h"
#include "positions.h"
#include "rand_pos.h"

using namespace std;

int main(int argc, char* argv[]) {
    //(void)argc;
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <input_name> <output_name> <rc> <kT> <n_insert> <seed>" << endl;
        return 1;  // Exit with error code indicating incorrect usage
    }

    string output_name = argv[1];
    int total_n_atoms = stoi(argv[3]);
    double T = stod(argv[2]);
    double rho = stod(argv[3]);
    int seed = stoi(argv[4]);
    int num_trails = stoi(argv[5]);

    // Start time
    // auto start_time = chrono::high_resolution_clock::now();
    // auto start_time_str = chrono::system_clock::to_time_t(start_time);
    // cout << "Start time: " << put_time(localtime(&start_time_str), "%Y-%m-%d %X") << endl;

    // Declare variables
    vector<vector<double> > box_dim(3, vector<double>(3));
    box_dim[0][0] = 15;
    box_dim[0][1] = 0;
    box_dim[0][2] = 0;
    box_dim[1][0] = 0;
    box_dim[1][1] = 15;
    box_dim[1][2] = 0;
    box_dim[2][0] = 0;
    box_dim[2][1] = 0;
    box_dim[2][2] = 15;  
    
    vector<vector<double> > positions(total_n_atoms, vector<double>(3));
    
    const double k = 8.617e-5;
    const double beta = 1/(k*T);
    const double e = (k*T)/1.2;
    const double V = box_dim[0][0]*box_dim[1][1]*box_dim[2][2];
    const double s_3 = rho*V/(total_n_atoms-1);
    const double s = pow(s_3, 1.0/3.0);

    // generate positions
    generateParticles(positions, total_n_atoms, rho, box_dim);

    //within the loop
    int trials = 0;
    int total_trials = 0;
    double step_size = 1.0;
    int accepted_moves = 0;
    int total_accepted_moves = 0;
    double w = 0;
    //till the insertion happens
    while(total_accepted_moves < num_trails){
        
        int i = rand() % total_n_atoms;

        // Calculate the initial energy of atom i
        double en_0 = e_i(0, e, box_dim, s, positions, i, total_n_atoms);

        double x_old = positions[i][0];
        double y_old = positions[i][1];
        double z_old = positions[i][2];

        // Generate a new position for atom i
        Coordinates new_pos = pos(box_dim, step_size);
        double x_new = new_pos.x;
        double y_new = new_pos.y;
        double z_new = new_pos.z;        
        
        positions[i][0] = x_new;
        positions[i][1] = y_new;
        positions[i][2] = z_new;

        // Calculate the new energy of atom i
        double en_new = e_i(0, e, box_dim, s, positions, i, total_n_atoms);

        uniform_real_distribution<> dis_real(0.0, 1.0);
        mt19937 gen(seed);
        double R = dis_real(gen);

        // Metropolis acceptance criterion
        if (R < exp(-beta * (en_new - en_0))) {
            // Accept the move, position is already updated
            w += exp(-beta * en_new);
            accepted_moves++;//register the insertion
            total_accepted_moves++;
        } 
        positions[i][0] = x_old;
        positions[i][1] = y_old;
        positions[i][2] = z_old;
            
        trials++;
        total_trials++;
        
        double acceptance_ratio = static_cast<double>(accepted_moves) / trials;
        cout << "Step: " << total_trials << ", Acceptance Ratio: " << acceptance_ratio << endl;

        if (acceptance_ratio < 0.3) {
            step_size *= 0.9; // Decrease step size
            // Reset counters
            trials = 0;
            accepted_moves = 0;
        } else if (acceptance_ratio > 0.5) {
            step_size *= 1.1; // Increase step size
            // Reset counters
            trials = 0;
            accepted_moves = 0;
        }    
    }
    
    //cout << trials << " number of trials " << endl; 
    cout << total_trials << " total number of trials " << endl;
    double acceptance_ratio = static_cast<double>(total_accepted_moves) / total_trials;
    cout << "avg Acceptance Ratio: " << acceptance_ratio << endl;
    cout << "w" << w << endl;
    
    // End time
    // auto end_time = chrono::high_resolution_clock::now();
    // auto end_time_str = chrono::system_clock::to_time_t(end_time);
    // cout << "End time: " << put_time(localtime(&end_time_str), "%Y-%m-%d %X") << endl;

    // Print processing time
    // chrono::duration<double> elapsed_time = end_time - start_time;
    // cout << "Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;

    return 0;
}
//distances need to be scaled
//PE has to be scaled
//density has to be scaled
//set the isotherm
//set the density


//analytical vs the code PE
//equation of state vs the test particle insertion method

