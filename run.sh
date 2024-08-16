#!/bin/bash

# Compile the program
g++ -std=c++11 -g -o mc_main mc_main.cpp pe_i.cpp generate_positions.cpp random.cpp output_func.cpp mc_eq.cpp pe_total.cpp mc_move.cpp

# Define the base input values
output_name="CONTCAR"
total_n_atoms=100
T=4.0
rho=0.01
seed=111
num_moves=200

csv_file="values.csv"

echo "ρ*,P*,μ_ex,acc ratio,Eq time (s), Wp time (S)"

# Loop to run the program multiple times with incrementally different input values
for i in {1..70}; do
    # Modify the input values incrementally
    current_rho=$(echo "$rho + 0.01 * $i" | bc)

    # Run the program with the modified input values
    program_output=$(./mc_main "$output_name" "$total_n_atoms" "$T" "$current_rho" "$seed" "$num_moves")

    # Check if the program ran successfully
    if [ $? -eq 0 ]; then
        echo "$i,$output_name,$total_n_atoms,$T,$current_rho,$seed,$num_moves,\"$program_output\"" >> "$csv_file"
    else
        echo "Error: ./mc_main failed with inputs: $output_name $total_n_atoms $T $current_rho $seed $num_moves"
        exit 1
    fi

    # Optionally, you can add a delay between runs
    # sleep 1
done
