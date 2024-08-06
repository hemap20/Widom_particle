import numpy as np
import random
import math

# Set simulation parameters
num_particles = N
epsilon = 1.0  # Lennard-Jones epsilon
sigma = 1.0    # Lennard-Jones sigma
temperature = 1.0  # Desired temperature
num_steps = 10000  # Number of MC steps for equilibration
box_size = [10.0, 10.0, 10.0]  # Dimensions of the simulation box
cutoff_radius = 2.5 * sigma

# Initialize positions
positions = initialize_positions(num_particles, box_size)

# Function to calculate the Lennard-Jones potential
def lennard_jones_potential(r2, epsilon, sigma):
    r6 = (sigma ** 2 / r2) ** 3
    r12 = r6 * r6
    return 4 * epsilon * (r12 - r6)

# Function to calculate total potential energy
def total_energy(positions, epsilon, sigma, cutoff_radius):
    energy = 0.0
    num_particles = len(positions)
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            dx = positions[i][0] - positions[j][0]
            dy = positions[i][1] - positions[j][1]
            dz = positions[i][2] - positions[j][2]
            # Apply minimum image convention for periodic boundary conditions
            dx -= box_size[0] * round(dx / box_size[0])
            dy -= box_size[1] * round(dy / box_size[1])
            dz -= box_size[2] * round(dz / box_size[2])
            r2 = dx * dx + dy * dy + dz * dz
            if r2 < cutoff_radius ** 2:
                energy += lennard_jones_potential(r2, epsilon, sigma)
    return energy

# Function to perform a single Monte Carlo step
def monte_carlo_step(positions, epsilon, sigma, temperature, box_size, cutoff_radius):
    num_particles = len(positions)
    k_B = 1.0  # Boltzmann constant in reduced units
    # Randomly select a particle and propose a displacement
    i = random.randint(0, num_particles - 1)
    displacement = [random.uniform(-0.1, 0.1) for _ in range(3)]
    new_position = [(positions[i][d] + displacement[d]) % box_size[d] for d in range(3)]
    # Calculate energy difference
    old_energy = 0.0
    new_energy = 0.0
    for j in range(num_particles):
        if j != i:
            dx_old = positions[i][0] - positions[j][0]
            dy_old = positions[i][1] - positions[j][1]
            dz_old = positions[i][2] - positions[j][2]
            dx_new = new_position[0] - positions[j][0]
            dy_new = new_position[1] - positions[j][1]
            dz_new = new_position[2] - positions[j][2]
            dx_old -= box_size[0] * round(dx_old / box_size[0])
            dy_old -= box_size[1] * round(dy_old / box_size[1])
            dz_old -= box_size[2] * round(dz_old / box_size[2])
            dx_new -= box_size[0] * round(dx_new / box_size[0])
            dy_new -= box_size[1] * round(dy_new / box_size[1])
            dz_new -= box_size[2] * round(dz_new / box_size[2])
            r2_old = dx_old * dx_old + dy_old * dy_old + dz_old * dz_old
            r2_new = dx_new * dx_new + dy_new * dy_new + dz_new * dz_new
            if r2_old < cutoff_radius ** 2:
                old_energy += lennard_jones_potential(r2_old, epsilon, sigma)
            if r2_new < cutoff_radius ** 2:
                new_energy += lennard_jones_potential(r2_new, epsilon, sigma)
    delta_energy = new_energy - old_energy
    # Apply Metropolis criterion
    if delta_energy < 0.0 or random.uniform(0.0, 1.0) < math.exp(-delta_energy / (k_B * temperature)):
        positions[i] = new_position  # Accept the move
        return True
    else:
        return False  # Reject the move

# Equilibration phase
for step in range(num_steps):
    monte_carlo_step(positions, epsilon, sigma, temperature, box_size, cutoff_radius)
    # Monitor properties (e.g., total energy) to check for equilibration
    if step % 100 == 0:
        energy = total_energy(positions, epsilon, sigma, cutoff_radius)
        print(f"Step {step}, Total Energy: {energy}")

# After equilibration, collect data for analysis
for step in range(num_production_steps):
    monte_carlo_step(positions, epsilon, sigma, temperature, box_size, cutoff_radius)
    # Collect data for analysis
    collect_data(positions)
