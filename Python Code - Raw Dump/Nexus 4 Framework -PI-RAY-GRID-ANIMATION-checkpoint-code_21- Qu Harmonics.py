import numpy as np
import matplotlib.pyplot as plt

# Initialize the grid
grid_size = 10
grid = np.zeros((grid_size, grid_size))
initial_active = int(0.1 * grid_size**2)  # 10% active cells
active_indices = np.random.choice(grid_size**2, initial_active, replace=False)
grid.flat[active_indices] = 1

# Function to get neighbors
def get_neighbors(i, j, grid_size):
    neighbors = []
    if i > 0: neighbors.append((i-1, j))
    if i < grid_size-1: neighbors.append((i+1, j))
    if j > 0: neighbors.append((i, j-1))
    if j < grid_size-1: neighbors.append((i, j+1))
    return neighbors

# Simulation parameters
threshold = 0.35
steps = 10

# Run the simulation
for step in range(steps):
    new_grid = grid.copy()
    for i in range(grid_size):
        for j in range(grid_size):
            neighbors = get_neighbors(i, j, grid_size)
            num_neighbors = len(neighbors)
            active_neighbors = sum(grid[n] for n in neighbors)
            proportion_active = active_neighbors / num_neighbors if num_neighbors > 0 else 0
            if proportion_active > threshold:
                new_grid[i, j] = 1
            else:
                new_grid[i, j] = 0  # Cells revert to 0 if threshold not met
    grid = new_grid

# Display results
print("Initial state had 10 active cells. Final state:")
print(f"Active cells: {int(np.sum(grid))}")
plt.imshow(grid, cmap='binary')
plt.title('Final Universe State')
plt.show()