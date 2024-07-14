#ifndef OUTPUT_FUNC_H
#define OUTPUT_FUNC_H

#include <string>
#include <vector>

using namespace std;

void print_CONTCAR(const string& filename, const vector<string>& atom_name, int n_atom_types, int total_n_atoms,
                double value, const vector<vector<double>>& box_dim,const string& coordinate_sys, const vector<vector<double>>& positions);

#endif // OUTPUT_FUNCT_H
