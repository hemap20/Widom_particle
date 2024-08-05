#include "output_func.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

void print_CONTCAR( const string filename, int N, const vector<vector<double> > box_dim, const vector<vector<double> > positions) {
    
    ofstream output(filename);

    if (!output.is_open()) {
        cerr << "Error opening file for writing: " << filename << endl;
        return;
    }

    // Print box dimensions
    for (int i = 0; i < 3; ++i) {
        output << "\t" << box_dim[i][0] << " " << box_dim[i][1] << " " << box_dim[i][2] << endl;
    }
    // Print positions of atoms
    for (int i = 0; i < N; ++i) {
        output << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;
    }
    output << endl;

    output.close();
}
