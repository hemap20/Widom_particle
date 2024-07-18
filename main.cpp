#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include <random>
#include "input_func.h"
#include "output_func.h"
#include "pairwise_dist.h"
#include "pot_energy.h"
#include "forces.h"
#include "insert.h"

using namespace std;

int main(int argc, char* argv[]) {
    //(void)argc;
    if (argc < 7) {
        cerr << "Usage: " << argv[0] << " <input_name> <output_name> <rc> <kT> <n_insert> <seed>" << endl;
        return 1;  // Exit with error code indicating incorrect usage
    }
    string input_name = argv[1];
    string output_name = argv[2];
    double rc = stod(argv[3]);
    double kT = stod(argv[4]);
    int n_insert = stoi(argv[5]);
    int seed = stoi(argv[6]);
    //bool print_flag = (stoi(argv[7]) != 0);

    // // Start time
    // auto start_time = chrono::high_resolution_clock::now();
    // auto start_time_str = chrono::system_clock::to_time_t(start_time);
    // cout << "Start time: " << put_time(localtime(&start_time_str), "%Y-%m-%d %X") << endl;

    // Declare variables
    vector<string> atom_name;
    int n_atom_types = 0;
    int total_n_atoms = 0;
    double value = 0;
    vector<vector<double> > box_dim(3, vector<double>(3));
    vector<int> n_atoms_per_type;
    string coordinate_sys;
    vector<vector<double> > positions;
    vector<double> distances;
    vector<tuple<int, int, double, vector<PairwiseDistance> > > pairwise_distances;
    vector<tuple<int, double, vector<PairwiseForce> > > pairwise_forces;    
    double beta = 1/kT;

    //generate a random number
    mt19937 gen(seed);

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    //compute distances
    dist(total_n_atoms, rc, box_dim, positions, pairwise_distances);

    //PE for the current configuration
    double PE_old = 0;
    PE_old = pot_energy(pairwise_distances, rc);
    cout << "PE_old " << PE_old << endl;

    //within the loop
    int n_acc = 0;
    int trials = 5;
    //till the insertion happens
    //while(n_acc<n_insert){
    for(int i=0; i<trials; i++){
        
        //perform insertion
        cout << "inserting main.cpp" << endl;
        insert_atom(total_n_atoms, box_dim, positions);

        //compute updated distances
        cout << "updated distances main.cpp" << endl;
        pairwise_distances.clear();
        dist(total_n_atoms, rc, box_dim, positions, pairwise_distances);

        //PE for current configuration
        double PE_new = 0;while(n_acc<n_insert)
        PE_new = pot_energy(pairwise_distances, rc);
    
        uniform_real_distribution<> dis_real(0.0, 1.0);

        //conditionally accept
        if(dis_real(gen) < exp(-beta*(PE_new-PE_old))){ 
            n_acc++; //register the insertion
            cout<< "PE_new " << PE_new << endl;
            cout << "registered" << endl;
            PE_old = PE_new;
        }
        else{
            //revert to the original positions, pairwise dist, total_num
            if (!positions.empty()) {
                cout << "popping back main.cpp" << endl;
                positions.pop_back();
                total_n_atoms = positions.size();
            }
            cout << "back to old dist main.cpp" << endl;
        }
        //trials++;
    }

    /* outside the loop: PE_old: find the current energy of the system
    inside the loop: insert the atom
    PE_new: find the energy of the new system
    check if the energy diff isn't too high, PE_old = PE_new
    else repeat the insertion, revert the changes to the dist, pos data strs
    */
    
    cout << trials << " number of trials " << endl; 
    //print the updated contcar
    print_CONTCAR(output_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, coordinate_sys, positions);
    
    // End time
    // auto end_time = chrono::high_resolution_clock::now();
    // auto end_time_str = chrono::system_clock::to_time_t(end_time);
    // cout << "End time: " << put_time(localtime(&end_time_str), "%Y-%m-%d %X") << endl;

    // Print processing time
    // chrono::duration<double> elapsed_time = end_time - start_time;
    // cout << "Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;

    return 0;
}


