#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include <random>
#include "pe_i.h"
#include "generate_positions.h"
#include "random.h"
#include "output_func.h"
#include "mc_eq.h"
#include "pe_total.h"

using namespace std;

int main(int argc, char* argv[]) {
    //(void)argc;
    if (argc < 6) {
        cerr << "Usage: " << argv[0] << " <output_name> <total_n_atoms> <T> <rho> <seed>" << endl;
        return 1;  // Exit with error code indicating incorrect usage
    }

    string output_name = argv[1];
    int total_n_atoms = stoi(argv[2]);
    double T = stod(argv[3]);
    double rho = stod(argv[4]);
    int seed = stoi(argv[5]);
    

    // Declare variables
    vector<vector<double> > box_dim(3, vector<double>(3));
    box_dim[0][0] = 20;
    box_dim[0][1] = 0;
    box_dim[0][2] = 0;
    box_dim[1][0] = 0;
    box_dim[1][1] = 20;
    box_dim[1][2] = 0;
    box_dim[2][0] = 0;
    box_dim[2][1] = 0;
    box_dim[2][2] = 20;  
    
    vector<vector<double> > positions(total_n_atoms, vector<double>(3));
    
    const double k = 8.617e-5; // Boltzmann constant in eV/K
    const double beta = 1/(k*T); // Inverse temperature
    const double e = (k*T)/1.2; // epsilon in eV
    const double V = box_dim[0][0]*box_dim[1][1]*box_dim[2][2]; 
    const double s_3 = rho*V/(total_n_atoms-1); 
    const double s = pow(s_3, 1.0/3.0); //sigma in Angstrom

    // generate positions
    generate_positions(positions, total_n_atoms, rho, box_dim);

    //initial energy
    double E = total_e(e, box_dim, s, positions, total_n_atoms);

    //print the initial positions
    print_CONTCAR(output_name, total_n_atoms, box_dim, positions);

    //within the loop
    int trials = 0;
    double step_size = 1.0;
    int accepted_moves = 0;
    double w = 0;
    
    ofstream PE("PE.csv");
    PE << "PE, time" << endl;

    // Record the start time
    // Start time
    auto start_time = chrono::system_clock::now();
    auto start_time_str = chrono::system_clock::to_time_t(start_time);

    //equilibration
    int total_trials = 10000;
    int total_accepted_moves = 0;
    for(int j = 0; j < total_trials; j++) {
        mc_eq(e, beta, s, box_dim, positions, total_n_atoms, w, step_size, trials, accepted_moves, total_accepted_moves, seed, E);  
        // Record the current time and calculate elapsed time
        auto current_time = chrono::system_clock::now();
        chrono::duration<double> elapsed = current_time - start_time;
        PE << E << "," << elapsed.count() << endl;
    }
    
    print_CONTCAR("EQUBM", total_n_atoms, box_dim, positions);

    //cout << trials << " number of trials " << endl; 
    // cout << total_trials << " total number of trials " << endl;
    // double acceptance_ratio = static_cast<double>(total_accepted_moves) / total_trials;
    // cout << "avg Acceptance Ratio: " << acceptance_ratio << endl;
    cout << "avg w: " << w/total_trials << endl;
    
    //print start time
    cout << "Start time: " << put_time(localtime(&start_time_str), "%Y-%m-%d %X") << endl;

    // End time
    auto end_time = chrono::system_clock::now();
    auto end_time_str = chrono::system_clock::to_time_t(end_time);
    cout << "End time: " << put_time(localtime(&end_time_str), "%Y-%m-%d %X") << endl;
    
    // Print processing time
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout << "Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;


    return 0;
}

//distances need to be scaled
//PE has to be scaled
//density has to be scaled
//set the isotherm
//set the density


//analytical vs the code PE
//equation of state vs the test particle insertion method
//the pressure of the system must change to accomodate the inserted atom
//new position is made, not an insertion

