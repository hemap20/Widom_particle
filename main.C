#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include "input_func.h"
#include "output_func.h"
#include "find_dist.h"

using namespace std;

int main(int argc, char* argv[]) {

    string input_name = argv[1];
    string output_name = argv[2];
    int ith_particle = stoi(argv[3]);

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

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    //compute distances
    find_dist(total_n_atoms, box_dim, positions, ith_particle, distances);

    // Print CONTCAR
    print_CONTCAR(output_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions, distances);

    // End time
    auto end_time = chrono::high_resolution_clock::now();
    auto end_time_str = chrono::system_clock::to_time_t(end_time);
    cout << "End time: " << put_time(localtime(&end_time_str), "%Y-%m-%d %X") << endl;

    // Print processing time
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout << "Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;

    return 0;
}

//g++ -o main main.C input_funct.C output_func.C -std=c++11
//./main input_file output_file ith_particle


