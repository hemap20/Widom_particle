#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <vector>
#include <string>

int main(){
    
    std::string line;
    int n_atom_types=0;
    double box_dim[3][3] = {0}; //Lx, Ly, Lz

    
    std::ifstream input("POSCAR");

    std::getline(input, line);   //KCl
    std::istringstream atoms(line);

    std::string token; //tokens of line stored here
    std::vector<std::string> atom_name = {};  
    
    //the number of types of atoms present
    while (atoms >> token) {
        n_atom_types++;
        atom_name.push_back(token);
    }
    token.clear();


    //1.0
    std::getline(input, line);   
    std::istringstream val(line);
    double value = 0;
    val >> value;


    //box dimensions
    for(int i=0; i<3; i++){
        std::getline(input, line);
        std::istringstream box(line);
        box >> box_dim[i][0] >> box_dim[i][1] >> box_dim[i][2];
    }

    //number of atoms of each type
    std::getline(input, line);
    std::istringstream atom_type(line);
    int n_atoms_per_type[n_atom_types];
    int total_n_atoms = 0;
    for(int i=0; i<n_atom_types; i++){
        int token_value;
        atom_type >> token_value;
        n_atoms_per_type[i] = token_value;
        total_n_atoms += token_value;
    }

    //coordinate system   
    std::string coordinate_sys;
    std::getline(input, coordinate_sys);

    //coordinates of the particles
    double positions[total_n_atoms][3];
    for(int i=0; i<total_n_atoms; i++){
        std::getline(input, line);
        std::istringstream coord(line);
        coord >> positions[i][0] >> positions[i][1] >> positions[i][2];
    }

    //closing input
    input.close();

    //writing CONTCAR file
    std::ofstream output("CONTCAR");

    //atom names
    for(int i=0; i<n_atom_types; i++){
       output <<  atom_name[i] << " ";
    }
    output << " " << std::endl;

    //charges
    output << value << std::endl;

    //box dimensions
    for(int i=0; i<3; i++){
        output << "\t" << box_dim[i][0] << " " << box_dim[i][1] << " " << box_dim[i][2] << std::endl;
    }

    //number of atoms per type
    for(int i=0; i<n_atom_types; i++){
        output << n_atoms_per_type[i] << " ";
    }
    output << std::endl;

    //coordinate system
    output << coordinate_sys << std::endl;

    //positions of atoms
    for(int i=0; i<total_n_atoms; i++){
        output << positions[i][0] << " " << positions[i][1] << " " << positions[i][2] << std::endl;
    }

    output.close();

    return 0;
}
