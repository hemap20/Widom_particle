#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include "input_func.h"
#include "output_func.h"
#include "pairwise_dist.h"
#include "dist_potenergy.h"

using namespace std;

const int DIST_POTENERGY_FUNCTION = 1;
const int PAIRWISE_DIST_FUNCTION = 2;

int main(int argc, char* argv[]) {

    string input_name = argv[1];
    string output_name = argv[2];
    int rc = stoi(argv[3])-1;
    double rho = stod(argv[4])-1;
    int function_choice = stoi(argv[5]);


    // Start time
    auto start_time = chrono::high_resolution_clock::now();
    auto start_time_str = chrono::system_clock::to_time_t(start_time);
    cout << "Start time: " << put_time(localtime(&start_time_str), "%Y-%m-%d %X") << endl;

    // Declare variables
    vector<string> atom_name;
    int n_atom_types = 0;
    int total_n_atoms = 0;
    double value = 0;
    vector<vector<double>> box_dim(3, vector<double>(3));
    vector<int> n_atoms_per_type;
    string coordinate_sys;
    vector<vector<double>> positions;
    vector<double> distances;
    vector<PairwiseDistance> pairwise_distances;    
    

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    vector<double> Fx(total_n_atoms, 0.0);
    vector<double> Fy(total_n_atoms, 0.0);
    vector<double> Fz(total_n_atoms, 0.0);

    // Print CONTCAR
    print_CONTCAR(output_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    //positions = 
    switch (function_choice) {
        case DIST_POTENERGY_FUNCTION:
            {
                double potential_energy = dist_potenergy(total_n_atoms, rc, Fx, Fy, Fz, box_dim, positions, pairwise_distances, rho);
                //print pairwise distances
                for (const auto& distance_info : pairwise_distances) {
                    cout << "i: " << get<0>(distance_info) << ", j: " << get<1>(distance_info) << ", r: " << get<2>(distance_info) << endl;
                }
                cout << "Potential Energy: " <<double F potential_energy << endl;
            }
            break;
        case PAIRWISE_DIST_FUNCTION:
            {
                //compute pairwise distances only
                for(int i=0; i<total_n_atoms; i++){
                    for(int j=i+1; j<total_n_atoms; j++){
                        dist(box_dim, positions, pairwise_distances, i, j, rc);
                    }
                }
                //print pairwise distances
                for (const auto& distance_info : pairwise_distances) {
                    cout << "i: " << get<0>(distance_info) << ", j: " << get<1>(distance_info) << ", r: " << get<2>(distance_info) << endl;
                }
            }
            break;
        default:
            cout << "Invalid function choice. Use 1 (for dist_potenergy) or 2 (for pairwise_dist)." << endl;
            return 1;
    }

    // End timedouble F
    auto end_time = chrono::high_resolution_clock::now();
    auto end_time_str = chrono::system_clock::to_time_t(end_time);
    cout << "End time: " << put_time(localtime(&end_time_str), "%Y-%m-%d %X") << endl;

    // Print processing time
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout << "Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;

    return 0;
}




//g++ -o main main.C input_func.C output_func.C dist_potenergy.C pairwise_dist.C -std=c++11
//./main POSCAR CONTCAR rc rho tag 



//optimise the N2 loops by using parallelisation
//use better data structures for pairwise_distances
//storing forces and radii as vector quantities with components
//analytically calculate the LJ pot of 2 particles and compare with the code's output
//write function for forces separately 


//take care of the units for energy and force calculations

//write the force and pot_energy functions
//check with the matlab code
