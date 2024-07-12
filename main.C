#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include <gsl/gsl_rng.h>
#include "input_func.h"
#include "output_func.h"
#include "pairwise_dist.h"
#include "pot_energy.h"
#include "forces.h"
#include "insert.h"

using namespace std;

int main(int argc, char* argv[]) {

    string input_name = argv[1];
    string output_name = argv[2];
    double rc = stod(argv[3]);
    double kT = stod(argv[4]);
    int n_insert = stoi(argv[5]);
    int seed = stoi(argv[6]);
    bool print_flag = (stoi(argv[7]) != 0);

    // Start time
    auto start_time = chrono::high_resolution_clock::now();
    auto start_time_str = chrono::system_clock::to_time_t(start_time);
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
    vector<double> distances;
    vector<tuple<int, int, double, vector<PairwiseDistance>>> pairwise_distances;
    vector<tuple<int, double, vector<PairwiseForce>>> pairwise_forces;    
    double beta = 1/kT;

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    //compute distances
    dist(total_n_atoms, rc, box_dim, positions, pairwise_distances);

    /*call pot energy
        find the potential energy for the current set of atoms    
    */
    double PE_old = pot_energy(pairwise_distances, rc);
    cout<< "PE_old" << PE_old << endl;


    //generate random r value
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r,seed);

    //pass total_num by reference

    /* insert particle
        generate a random set of coordinates from r
        add to the positions list, update total_number_atoms
    */
   insert_atom(r, total_n_atoms, box_dim, positions);

   //compute updated distances
    pairwise_distances.clear();
    dist(total_n_atoms, rc, box_dim, positions, pairwise_distances);

   /*call pot energy
        find the potential of the updated positions list
   */
    double PE_new = pot_energy(pairwise_distances, rc);
    cout<< "PE_new" << PE_new << endl;

  /*check for a valid insertion
        using the metropolis algo
        if valid: register the move; N++
        else: revert back to the original list
  */
    //within the loop
    // int n_acc = 0;
    // //till the insertion happens
    // while(n_acc < n_insert){
    //     //gsl_rng_uniform(r) has to be b/n 0,1
    //     if(gsl_rng_uniform(r) < exp(-beta*(PE_new-PE_old))){ //should depend on the density
    //         PE_old = PE_new;
    //         n_acc++; //register the insertion 
    //     }
    //     else{
    //         //revert to the original positions, pairwise dist, total_num
    //         positions.pop_back();
    //         //*********clear dist before this
    //         dist(total_n_atoms, rc, box_dim, positions, pairwise_distances);
    //         total_n_atoms = positions.size();
    //     }
    // }
    
    //print the updated contcar
    
    // End time
    auto end_time = chrono::high_resolution_clock::now();
    auto end_time_str = chrono::system_clock::to_time_t(end_time);
    cout << "End time: " << put_time(localtime(&end_time_str), "%Y-%m-%d %X") << endl;

    // Print processing time
    chrono::duration<double> elapsed_time = end_time - start_time;
    cout << "Processing time: " << fixed << setprecision(6) << elapsed_time.count() << " seconds" << endl;

    return 0;
}

//optimise the N2 loops by using parallelisation
//take an avegage in the end

