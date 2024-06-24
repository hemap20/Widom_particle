#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

int main(){
    
    std::string line;
    int n_atom_types=0;
    double box_dim[3] = {0}; //Lx, Ly, Lz

    
    std::ifstream input("POSCAR");

    getline(input, line);   //KCl
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
    getline(input, line);   

    //box dimensions
    for(int i=0; i<3; i++){
        int box_token = 0;
        getline(input, line);
        std::istringstream box(line);
        if(box >> box_token){
            box_dim[i] = box_token;
        }
    }

    //number of atoms of each type
    getline(input, line);
    std::istringstream atom_type(line);
    int n_atoms_per_type[n_atom_types] = {0};
    for(int i=0; i<n_atom_types; i++){
        
        
    }

     
    
    






    
}