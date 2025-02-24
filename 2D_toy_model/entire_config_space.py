import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Define the complex potential energy surface as a sum of several 2D Gaussians
def potential(x, y):
    gaussians = [
        {"center": (0.5, 0.5), "amplitude": -5, "width": 0.25},
        {"center": (-0.5, -0.5), "amplitude": -5, "width": 0.25},
        {"center": (0.0, 0.0), "amplitude": -5, "width": 0.25},
        {"center": (0.5, -0.5), "amplitude": -5, "width": 0.25},
        {"center": (-0.5, 0.5), "amplitude": -5, "width": 0.25},
    ]
    potential_value = 0
    for g in gaussians:
        xc, yc = g["center"]
        amplitude = g["amplitude"]
        width = g["width"]
        potential_value += amplitude * np.exp(-((x - xc)**2 + (y - yc)**2) / (2 * width**2))
    return potential_value

def degree_of_overlap(x,y):
    doo = 1
    if (-0.5 <= x < 0.5) and (-0.5 <= y < 0.5):
        doo = 2
    return doo

# Monte Carlo simulation parameters
n_steps = 10000000
step_size = 0.1
boundary = 1.0
boundary_xmin, boundary_xmax = -1, 1
boundary_ymin, boundary_ymax = -1, 1
x, y = np.random.uniform(boundary_xmin, boundary_xmax), np.random.uniform(boundary_ymin, boundary_ymax)  # Starting position
print("initial coordinate: (%.3f, %.3f)"%(x,y))

positions = [(x, y)]
potential_energies = [potential(x,y)]
DOOs = [degree_of_overlap(x,y)]

# Monte Carlo sampling loop
for _ in range(n_steps):
    # Propose a new position
    x_new = x + step_size * (np.random.rand() - 0.5)
    y_new = y + step_size * (np.random.rand() - 0.5)

    # Check boundary condition
    if x_new > boundary_xmax or x_new < boundary_xmin or y_new > boundary_ymax or y_new < boundary_ymin:
        continue

    # Metropolis criterion
    V_current = potential(x, y)
    V_new = potential(x_new, y_new)
    delta_potential =  V_new - V_current
    if delta_potential < 0 or np.random.rand() < np.exp(-delta_potential):
        x, y = x_new, y_new
        V_current = V_new

    positions.append((x, y))
    potential_energies.append(V_current)
    DOOs.append(degree_of_overlap(x,y))

# Extract x and y positions
positions = np.array(positions)
x_positions, y_positions = positions[:, 0], positions[:, 1]

# Visualize the potential energy surface and sampled positions
x_grid = np.linspace(-boundary, boundary, 200)
y_grid = np.linspace(-boundary, boundary, 200)
x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)
z_mesh = potential(x_mesh, y_mesh)

potential_energies = np.array(potential_energies)
DOOs = np.array(DOOs)
print("mean DOO = %.3f, mean potential energy = %.3f"%(np.mean(np.array(DOOs)), np.mean(np.array(potential_energies))))

""" plt.figure(figsize=(10, 8))
plt.contourf(x_mesh, y_mesh, z_mesh, levels=50, cmap="viridis")
plt.colorbar(label="Potential Energy")
plt.plot(x_positions, y_positions, 'r.', markersize=1, label="Sampled Positions")
plt.gca().set_aspect('equal', adjustable='box')
plt.axvline(boundary_xmin, color="black", linestyle="--")
plt.axvline(boundary_xmax, color="black", linestyle="--")
plt.axhline(boundary_ymin, color="black", linestyle="--")
plt.axhline(boundary_ymax, color="black", linestyle="--")
plt.xlim(-boundary - 0.05, boundary + 0.05)
plt.ylim(-boundary - 0.05, boundary + 0.05)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Monte Carlo Sampling on Complex 2D Potential")
plt.legend()
plt.show() """

