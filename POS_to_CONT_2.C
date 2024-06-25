#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>
#include <chrono> 
#include <sstream>

using namespace std;
using namespace chrono;

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

void print_CONTCAR(const string& filename, const vector<string>& atom_name, int n_atom_types, int total_n_atoms,
                double value, const vector<vector<double>>& box_dim, const vector<int>& n_atoms_per_type,
                const string& coordinate_sys, const vector<vector<double>>& positions) {
    
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

    // Print number of atoms per type
    for (int i = 0; i < n_atom_types; ++i) {
        output << n_atoms_per_type[i] << " ";
    }
    output << endl;

    // Print coordinate system
    output << coordinate_sys << endl;

    // Print positions of atoms
    for (int i = 0; i < total_n_atoms; ++i) {
        output << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;
    }

    output.close();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " input_file output_file" << endl;
        return 1;
    }

    string input_name = argv[1];
    string output_name = argv[2];

    //start time
    auto start_time = high_resolution_clock::now();
    auto start_time_str = system_clock::to_time_t(start_time);
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

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    // Print CONTCAR
    print_CONTCAR(output_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    //end time
    auto end_time = high_resolution_clock::now();
    auto end_time_str = system_clock::to_time_t(end_time);
    cout << "End time: " << put_time(localtime(&end_time_str), "%Y-%m-%d %X") << endl;

    //print time taken
    duration<double> elapsed_time = duration_cast<duration<double>>(end_time - start_time);
    cout << "Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;

    return 0;
}

//functions in different folders
//taking pairs of particles and computing the distance
