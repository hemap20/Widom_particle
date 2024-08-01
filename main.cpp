#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include <random>
#include "input_func.h"
#include "output_func.h"
#include "pairwise_dist.h"
#include "pot_energy.h"
#include "forces.h"
#include "insert.h"

using namespace std;

int main(int argc, char* argv[]) {
    //(void)argc;
    if (argc < 7) {
        cerr << "Usage: " << argv[0] << " <input_name> <output_name> <rc> <kT> <n_insert> <seed>" << endl;
        return 1;  // Exit with error code indicating incorrect usage
    }
    string input_name = argv[1];
    string output_name = argv[2];
    //double rc = stod(argv[3]);
    double T = stod(argv[3]);
    double rho = stod(argv[4]);
    int n_insert = stoi(argv[5]);
    int seed = stoi(argv[6]);

    // Start time
    // auto start_time = chrono::high_resolution_clock::now();
    // auto start_time_str = chrono::system_clock::to_time_t(start_time);
    // cout << "Start time: " << put_time(localtime(&start_time_str), "%Y-%m-%d %X") << endl;

    // Declare variables
    vector<string> atom_name;
    int n_atom_types = 0;
    int total_n_atoms = 0;
    double value = 0;
    vector<vector<double> > box_dim(3, vector<double>(3));
    vector<int> n_atoms_per_type;
    string coordinate_sys;
    vector<vector<double> > positions;
    vector<double> distances;
    vector<tuple<int, int, double, vector<PairwiseDistance> > > pairwise_distances;
    vector<tuple<int, double, vector<PairwiseForce> > > pairwise_forces;    
    
    const k=8.617e-5;
    double beta = 1/(k*T);
    const double e = (k*T)/1.2;

    //generate a random number
    mt19937 gen(seed);

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    const double V = box_dim[0][0]*box_dim[1][1]*box_dim[2][2];
    const double s_3 = rho*V/(total_n_atoms-1);
    const double s = pow(s_3, 1.0/3.0);

    //compute distances
    //dist(total_n_atoms, box_dim, positions, pairwise_distances);

    //PE for the current configuration
    double PE_old = 0;
    PE_old = total_e(e, s, total_n_atoms, box_dim, positions, pairwise_distances);
    //cout << "PE_old " << PE_old << endl;

    ofstream csvFile;
    csvFile.open("PE.csv");
    csvFile << "PE, step" << endl;
    //csvFile << PE_old << " ," << 0 <<endl;

    //within the loop
    int trials = 0;
    int total_trials = 0;
    double step_size = 1.0;
    int accepted_moves = 0;
    int total_accepted_moves = 0;
    //till the insertion happens
    for(int n_acc=0; n_acc<n_insert;){
        
        //perform insertion
        insert_atom(total_n_atoms, box_dim, positions, step_size);

        //compute updated distances
        // pairwise_distances.clear();
        // dist(total_n_atoms,box_dim, positions, pairwise_distances);

        //PE for current configuration
        double PE_new = 0;
        PE_new = e_i(total_n_atoms, s, total_n_atoms, box_dim, positions, pairwise_distances);
    
        uniform_real_distribution<> dis_real(0.0, 1.0);
        double R = dis_real(gen);
        //conditionally accept
        if( R/4 < exp(-beta*(PE_new-PE_old))){ 
            n_acc++; 
            accepted_moves++;//register the insertion
            total_accepted_moves++;
            //cout<< "PE_new " << PE_new << endl;
            csvFile << PE_new << " ," << total_trials << endl;
            PE_old = PE_new;
        }
        else{
            //revert to the original positions, pairwise dist, total_num
            if (!positions.empty()) {
                positions.pop_back();
                total_n_atoms = positions.size();
            }
            dist(total_n_atoms,box_dim, positions, pairwise_distances);
        }
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
    csvFile.close();
    
    //cout << trials << " number of trials " << endl; 
    cout << total_trials << " total number of trials " << endl;
    double acceptance_ratio = static_cast<double>(total_accepted_moves) / total_accepted_moves;
    cout << "avg Acceptance Ratio: " << acceptance_ratio << endl;
    //print the updated contcar
    print_CONTCAR(output_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, coordinate_sys, positions);
    
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

