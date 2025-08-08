import numpy as np

# Define the 3x3 grid
grid = np.zeros((8,8))
alpha_pos = (1, 1)  # Alpha at the center
grid[alpha_pos] = 1.0  # Alpha's initial strength

# Harmonic constant for decay
decay_constant = 0.35

# Calculate influence of alpha on each cell
for i in range(8):
    for j in range(8):
        if (i, j) != alpha_pos:
            distance = np.sqrt((i - alpha_pos[0])**2 + (j - alpha_pos[1])**2)
            grid[i, j] = np.exp(-decay_constant * distance)  # Influence decreases with distance

# Print the grid showing alpha's influence
print("Alpha's localized gravity influence on the 3x3 grid:")
print(grid)