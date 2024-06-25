#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip> 
#include <string>
#include <sstream>

using namespace std;

void read_input(const string& filename, vector<string>& atom_name, int& n_atom_types, int& total_n_atoms,
                double& value, double box_dim[3][3], int n_atoms_per_type[], string& coordinate_sys,
                double positions[][3]){
    
    ifstream input(filename);

    string line;

    //KCl
    getline(input, line);
    istringstream atoms(line);

    //tokens of line stored here 
    string token; 

    //the number of types of atoms present
    while (atoms >> token) {
        n_atom_types++;
        atom_name.push_back(token);
    }
    token.clear();


    //1.0
    getline(input, line);   
    istringstream val(line);
    val >> value;


    //box dimensions
    for(int i=0; i<3; i++){
        getline(input, line);
        istringstream box(line);
        box >> box_dim[i][0] >> box_dim[i][1] >> box_dim[i][2];
    }

    //number of atoms of each type
    getline(input, line);
    istringstream atom_type(line);
    for(int i=0; i<n_atom_types; i++){
        int token_value;
        atom_type >> token_value;
        n_atoms_per_type[i] = token_value;
        total_n_atoms += token_value;
    }

    //coordinate system   
    getline(input, coordinate_sys);

    //coordinates of the particles
    for(int i=0; i<total_n_atoms; i++){
        getline(input, line);
        istringstream coord(line);
        coord >> positions[i][0] >> positions[i][1] >> positions[i][2];
    }

    //closing input
    input.close();
}

void print_CONTCAR(const string& filename, vector<string>& atom_name, int& n_atom_types, int& total_n_atoms,
                double& value, double box_dim[3][3], int n_atoms_per_type[], string& coordinate_sys,
                double positions[][3]){
    
    ofstream output(filename);

    //atom names
    for(int i=0; i<n_atom_types; i++){
       output <<  atom_name[i] << " ";
    }
    output << " " << endl;

    //charges
    output << value << endl;

    //box dimensions
    for(int i=0; i<3; i++){
        output << "\t" << box_dim[i][0] << " " << box_dim[i][1] << " " << box_dim[i][2] << endl;
    }

    //number of atoms per type
    for(int i=0; i<n_atom_types; i++){
        output << n_atoms_per_type[i] << " ";
    }
    output << endl;

    //coordinate system
    output << coordinate_sys << endl;

    //positions of atoms
    for(int i=0; i<total_n_atoms; i++){
        output << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << endl;
    }

    output.close();
}


int main(int argc, char* argv[]){
    
    string input_name = argv[1];
    string output_name = argv[2];

    //all variables
    vector<string> atom_name = {};
    int n_atom_types=0;
    int total_n_atoms = 0;
    double value = 0;
    double box_dim[3][3] = {0};
    int n_atoms_per_type[n_atom_types];
    string coordinate_sys;
    double positions[total_n_atoms][3];

    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);
    print_CONTCAR(output_name, atom_name, n_atom_types,total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    return 0;
}