import numpy as np
import matplotlib.pyplot as plt
import math

# Define the function to calculate F_cubed with twist and spin contributions
def calculate_f_cubed_with_axes(A, B, len_c, x, axis_weights):
    twist = axis_weights[0] * (A**2)  # Twist driven by axis 1
    spin = axis_weights[1] * (B**2)  # Spin driven by axis 2
    torque = axis_weights[2] * (1 + math.exp(-10 * (A * x - 0.35)))  # Torque component
    return (twist + spin) * len_c * torque

# Recursive refinement function to calculate F_cubed with dynamic contributions
def recursive_f_cubed_with_axes(A, B, len_function, x, iterations):
    current_f = 0  # Initialize F_cubed value
    axis_weights = [1, 1, 1]  # Initial weights for axes
    values = []  # Store overall F_cubed values
    axis_contributions = []  # Store individual axis contributions

    for i in range(iterations):
        len_c = len_function(current_f)  # Calculate Len(C)
        axis_weights = [np.sin(i / 10), np.cos(i / 10), np.sin(i / 20)]  # Dynamic twist, spin, torque weights
        current_f = calculate_f_cubed_with_axes(A, B, len_c, x, axis_weights)  # Update F_cubed
        
        values.append(current_f)  # Store overall value
        axis_contributions.append(axis_weights)  # Store current axis weights

    return values, np.array(axis_contributions)

# Parameters for the calculation
A = 2.0
B = 3.0
x = 1.0
iterations = 2000

# Define dynamic Len(C) function for recursive refinement
len_function = lambda current_f: max(1, current_f + 1)  # Dynamic Len(C)

# Perform the recursive calculation
values, axis_contributions = recursive_f_cubed_with_axes(A, B, len_function, x, iterations)

# Plotting overall F_cubed values
plt.figure(figsize=(10, 6))
plt.plot(range(1, iterations + 1), values, label='F_cubed', color='b', marker='o', linestyle='-')
plt.title("Recursive Refinement of F_cubed Over Iterations", fontsize=14)
plt.xlabel("Iteration", fontsize=12)
plt.ylabel("F_cubed Value", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.show()

# Plotting contributions from each axis
plt.figure(figsize=(10, 6))
plt.plot(range(1, iterations + 1), axis_contributions[:, 0], label='Twist (Axis 1)', color='r', linestyle='--')
plt.plot(range(1, iterations + 1), axis_contributions[:, 1], label='Spin (Axis 2)', color='g', linestyle='--')
plt.plot(range(1, iterations + 1), axis_contributions[:, 2], label='Torque (Axis 3)', color='m', linestyle='--')
plt.title("Axis Contributions to Recursive Refinement", fontsize=14)
plt.xlabel("Iteration", fontsize=12)
plt.ylabel("Axis Weight (Dynamic Contribution)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.show()

# Energy conservation plot (sum of contributions over iterations)
energy_conservation = np.sum(axis_contributions, axis=1)
plt.figure(figsize=(10, 6))
plt.plot(range(1, iterations + 1), energy_conservation, label='Energy Conservation (Sum of Axes)', color='c')
plt.title("Energy Conservation Across Iterations", fontsize=14)
plt.xlabel("Iteration", fontsize=12)
plt.ylabel("Energy (Sum of Axes)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.show()
