#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <tuple>
#include "input_func.h"
#include "output_func.h"
#include "pairwise_dist.h"
#include "pot_energy.h"
#include "forces.h"

using namespace std;

int main(int argc, char* argv[]) {

    string input_name = argv[1];
    string output_name = argv[2];
    double rc = stod(argv[3])-1;
    double rho = stod(argv[4])-1;
    bool print_flag = (stoi(argv[5]) != 0);

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

    // Read input
    read_input(input_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    // Print CONTCAR
    print_CONTCAR(output_name, atom_name, n_atom_types, total_n_atoms, value, box_dim, n_atoms_per_type, coordinate_sys, positions);

    //compute distances
    dist(total_n_atoms, rc, box_dim, positions, pairwise_distances);

    //print the elements (optional)
    if(print_flag){
        for (const auto& item : pairwise_distances) {
            int i = get<0>(item);
            int j = get<1>(item);
            double r = get<2>(item);
            cout << "i = " << i << ", j = " << j << ", r = " << r << endl;
        }
    }

    //compute total potential energy
    cout << "Total Potential Energy: " << pot_energy(pairwise_distances, rc) << endl;

    //compute the forces
    forces(pairwise_distances, pairwise_forces);

    //print forces (optional)
    if(print_flag){
        for (const auto& item : pairwise_forces) {
            int i = get<0>(item);
            double F = get<1>(item);
            const vector<PairwiseForce>& forces = get<2>(item);
            cout << "Particle " << i << " Force Mag: " << F << endl;
            for (const auto& force : forces) {
                cout << "  Force vector: (";
                for (size_t k = 0; k < force.F_vec.size(); ++k) {
                    cout << force.F_vec[k];
                    if (k < force.F_vec.size() - 1) {
                        cout << ", ";
                    }
                }
                cout << ")\n";
            }
        }
    }

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
