#include "input_func.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

void read_input(const string& filename, vector<string>& atom_name, int& n_atom_types, int& total_n_atoms,
                double& value, vector<vector<double>>& box_dim, vector<int>& n_atoms_per_type,
                string& coordinate_sys, vector<vector<double>>& positions) {
    
    ifstream input(filename);

    if (!input.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    string line;

    // Read atom types
    getline(input, line);
    istringstream atoms(line);
    string token;
    while (atoms >> token) {
        n_atom_types++;
        atom_name.push_back(token);
    }

    // Read value
    getline(input, line);
    istringstream val(line);
    val >> value;

    // Read box dimensions
    box_dim.resize(3, vector<double>(3));
    for (int i = 0; i < 3; ++i) {
        getline(input, line);
        istringstream box(line);
        box >> box_dim[i][0] >> box_dim[i][1] >> box_dim[i][2];
    }

    // Read number of atoms per type
    getline(input, line);
    istringstream atom_type(line);
    int total_atoms = 0;
    for (int i = 0; i < n_atom_types; ++i) {
        int token_value;
        atom_type >> token_value;
        n_atoms_per_type.push_back(token_value);
        total_atoms += token_value;
    }
    total_n_atoms = total_atoms;

    // Read coordinate system
    getline(input, coordinate_sys);

    // Read positions
    positions.resize(total_n_atoms, vector<double>(3));
    for (int i = 0; i < total_n_atoms; ++i) {
        getline(input, line);
        istringstream coord(line);
        coord >> positions[i][0] >> positions[i][1] >> positions[i][2];
    }

    input.close();
}
