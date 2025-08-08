import numpy as np
import matplotlib.pyplot as plt

# Initialize fluid particles
num_particles = 100
positions = np.random.rand(num_particles, 2)  # X, Y positions
velocities = np.random.rand(num_particles, 2) - 0.5  # Random velocity

# Simulation parameters
num_steps = 500
positions_history = []

# Simulate particle motion
for step in range(num_steps):
    # Apply random noise to velocities (turbulence)
    turbulence = np.random.normal(scale=0.01, size=velocities.shape)
    velocities += turbulence

    # Update positions
    positions += velocities
    positions_history.append(positions.copy())

# Visualize fluid motion
plt.figure(figsize=(10, 6))
for snapshot in positions_history[::50]:  # Plot every 50 steps
    plt.scatter(snapshot[:, 0], snapshot[:, 1], s=10)
plt.title("Oscillatory Fluid Dynamics with Noise")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
