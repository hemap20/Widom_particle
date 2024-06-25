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

