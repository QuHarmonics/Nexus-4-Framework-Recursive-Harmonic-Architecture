import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# Define the parameters for the Navier-Stokes equations
rho = 1.0  # fluid density
nu = 0.01  # fluid viscosity
L = 1.0  # length of the domain

# Define the grid size and time step
N = 100
dt = 0.001

# Initialize the velocity and pressure fields
u = np.zeros((N, N))
v = np.zeros((N, N))
p = np.zeros((N, N))

# Define the boundary conditions
u[0, :] = 1.0
u[-1, :] = 0.0
v[0, :] = 0.0
v[-1, :] = 0.0

# Add a small perturbation to the initial velocity field
u[50, 50] = 0.1

# Time-stepping loop
for n in range(1000):
    # Compute the convective term
    conv_u = u * np.gradient(u)[0] / L
    conv_v = v * np.gradient(v)[0] / L

    # Compute the viscous term
    visc_u = nu * np.gradient(np.gradient(u)[0])[0] / L**2
    visc_v = nu * np.gradient(np.gradient(v)[0])[0] / L**2

    # Update the velocity field
    u += dt * (-conv_u + visc_u)
    v += dt * (-conv_v + visc_v)

    # Compute the pressure field
    p = np.linalg.solve(np.eye(N) - dt * nu * np.gradient(np.eye(N))[0], -rho * np.gradient(u)[0])

    # Apply a small amount of diffusion to the velocity field
    u += 0.01 * dt * np.gradient(np.gradient(u)[0])[0] / L**2
    v += 0.01 * dt * np.gradient(np.gradient(v)[0])[0] / L**2

# Plot the velocity and pressure fields
plt.figure(figsize=(12, 6))
plt.subplot(121)
plt.imshow(u, cmap='viridis')
plt.title('Velocity Field')
plt.subplot(122)
plt.imshow(p, cmap='viridis')
plt.title('Pressure Field')
plt.show()

# Plot the velocity field at different time steps
plt.figure(figsize=(12, 6))
for n in range(0, N, 10):  
    plt.plot(u[:, n], label=f'Grid point {n}')
plt.xlabel('Grid point')
plt.ylabel('Velocity')
plt.title('Velocity Field at Different Grid Points')
plt.legend()
plt.show()