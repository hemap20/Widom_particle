#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include "input_func.h"
#include "output_func.h"
#include "pairwise_dist.h"
//#include "dist_potenergy.h"
#include "pot_energy.h"
#include "forces.h"

using namespace std;

int main(int argc, char* argv[]) {

    string input_name = argv[1];
    string output_name = argv[2];
    int rc = stoi(argv[3])-1;
    double rho = stod(argv[4])-1;

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
    vector<tuple<int, int, double, vector<PairwiseDistance>>> pairwise_distances;
    vector<tuple<int,vector<PairwiseForce>>> pairwise_forces;    

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    vector<double> Fx(total_n_atoms, 0.0);
    vector<double> Fy(total_n_atoms, 0.0);
    vector<double> Fz(total_n_atoms, 0.0);

    // Print CONTCAR
    print_CONTCAR(output_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    //compute distances
    dist(total_n_atoms, rc, box_dim, positions, pairwise_distances, rho);

    //print the elements (optional)
    //***

    //compute total potential energy
    cout << "Total Potential Energy" << pot_energy(pairwise_distances, rc) << endl;

    //compute the forces
    forces(pairwise_distances, pairwise_forces);

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
//analytically calculate the LJ pot of 2 particles and compare with the code's output
//check with the matlab code
//change main to remove dist, include forces file
