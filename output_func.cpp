#include "output_func.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

void print_CONTCAR( const string& filename, const vector<string>& atom_name, int n_atom_types, int total_n_atoms,
                double value, const vector<vector<double>>& box_dim, const string& coordinate_sys, const vector<vector<double>>& positions) {
    
    ofstream output(filename);

    if (!output.is_open()) {
        cerr << "Error opening file for writing: " << filename << endl;
        return;
    }

    // Print atom names
    for (int i = 0; i < n_atom_types; ++i) {
       output << atom_name[i] << " ";
    }
    output << endl;

    // Print value
    output << value << endl;

    // Print box dimensions
    for (int i = 0; i < 3; ++i) {
        output << "\t" << box_dim[i][0] << " " << box_dim[i][1] << " " << box_dim[i][2] << endl;
    }

    // // Print number of atoms per type
    // for (int i = 0; i < n_atom_types; ++i) {
    //     output << n_atoms_per_type[i] << " ";
    // }
    // output << endl;

    // Print coordinate system
    output << coordinate_sys << endl;

    // Print positions of atoms
    total_n_atoms = positions.size();
    for (int i = 0; i < total_n_atoms+1; ++i) {
        output << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;
    }
    output << endl;

    output.close();
}
