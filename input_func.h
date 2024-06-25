#ifndef INPUT_FUNC_H
#define INPUT_FUNC_H

#include <string>
#include <vector>

using namespace std;

void read_input(const string& filename, vector<string>& atom_name, int& n_atom_types, int& total_n_atoms,
                double& value, vector<vector<double>>& box_dim, vector<int>& n_atoms_per_type,
                string& coordinate_sys, vector<vector<double>>& positions);

#endif // INPUT_FUNCT_H
