import numpy as np
from scipy.integrate import dblquad

# Define the 2D potential energy surface: summation of multiple Gaussians
def potential_energy(x, y):
    gaussian1 = -5 * np.exp(-((x - 0.5)**2 + (y - 0.5)**2) / (2 * 0.25**2))
    gaussian2 = -5 * np.exp(-((x + 0.5)**2 + (y + 0.5)**2) / (2 * 0.25**2))
    gaussian3 = -5 * np.exp(-(x**2 + y**2) / (2 * 0.25**2))
    gaussian4 = -5 * np.exp(-((x - 0.5)**2 + (y + 0.5)**2) / (2 * 0.25**2))
    gaussian5 = -5 * np.exp(-((x + 0.5)**2 + (y - 0.5)**2) / (2 * 0.25**2))
    return gaussian1 + gaussian2 + gaussian3 + gaussian4 + gaussian5

# Define the Boltzmann factor for the potential at temperature T
def boltzmann_factor(x, y, beta):
    return np.exp(-beta * potential_energy(x, y))

def potential_energy_boltzmann_factor(x, y, beta):
    return potential_energy(x,y) * np.exp(-beta * potential_energy(x, y))

def degree_of_overlap(x,y):
    doo = 1
    if (-0.5 <= x < 0.5) and (-0.5 <= y < 0.5):
        doo = 2
    return doo

def degree_of_overlap_boltzmann_factor(x, y, beta):
    return degree_of_overlap(x, y) * np.exp(-beta * potential_energy(x, y))

def potential_energy_tilde_boltzmann_factor(x, y, beta):
    return potential_energy(x, y) / degree_of_overlap(x, y) * np.exp(-beta * potential_energy(x, y))

def degree_of_overlap_reciprocal_boltzmann_factor(x, y, beta):
    return 1.0/degree_of_overlap(x, y) * np.exp(-beta * potential_energy(x, y))

# Define the integration region (rectangular boundary)
x_min, x_max = -1, 1
y_min, y_max = -1, 1

# Inverse temperature (beta = 1 / (k_B * T))
beta = 1.0  # Assuming k_B * T = 1 for simplicity

# Integrate the partition function Z
partition_function, error = dblquad(
    lambda y, x: boltzmann_factor(x, y, beta),  # Function to integrate
    x_min, x_max,  # Integration bounds for x
    lambda x: y_min,  # Lower bound for y (as a function of x)
    lambda x: y_max   # Upper bound for y (as a function of x)
)

average_potential_energy, error = dblquad(
    lambda y, x: potential_energy_boltzmann_factor(x, y, beta),  # Function to integrate
    x_min, x_max,  # Integration bounds for x
    lambda x: y_min,  # Lower bound for y (as a function of x)
    lambda x: y_max   # Upper bound for y (as a function of x)
)

print(f"average potential energy: {average_potential_energy/partition_function:.6f}")
print(f"Estimated Integration Error: {error:.6e}")

average_degree_of_overlap, error = dblquad(
    lambda y, x: degree_of_overlap_boltzmann_factor(x, y, beta),  # Function to integrate
    x_min, x_max,  # Integration bounds for x
    lambda x: y_min,  # Lower bound for y (as a function of x)
    lambda x: y_max   # Upper bound for y (as a function of x)
)

print(f"average degree of overlap: {average_degree_of_overlap/partition_function:.6f}")
print(f"Estimated Integration Error: {error:.6e}")
