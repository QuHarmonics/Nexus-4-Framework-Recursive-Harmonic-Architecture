import matplotlib.pyplot as plt
import math

# Define the function to calculate F_cubed
def calculate_f_cubed(A, B, len_c, x):
    # Quadratic component
    quadratic_part = A**2 + B**2
    
    # Exponential component
    exponential_part = 1 + math.exp(-10 * (A * x - 0.35))
    
    # Combine all components
    return quadratic_part * len_c * exponential_part

# Recursive refinement function for F_cubed
def recursive_f_cubed(A, B, len_function, x, iterations):
    current_f = 0  # Initialize current F_cubed value
    values = []  # To store F_cubed values at each iteration
    
    for _ in range(iterations):
        len_c = len_function(current_f)  # Dynamically calculate Len(C)
        current_f = calculate_f_cubed(A, B, len_c, x)  # Update F_cubed value
        values.append(current_f)  # Store the value for plotting
    
    return values

# Parameters for the calculation
A = 2.0
B = 3.0
x = 1.0
iterations = 200

# Define a dynamic Len(C) function
len_function = lambda current_f: max(1, current_f + 1)  # Example function

# Perform the recursive calculation
values = recursive_f_cubed(A, B, len_function, x, iterations)

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(range(1, iterations + 1), values, marker='o', linestyle='-', color='b', label='F_cubed')
plt.title("Recursive Refinement of F_cubed Over Iterations", fontsize=14)
plt.xlabel("Iteration", fontsize=12)
plt.ylabel("F_cubed Value", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=12)
plt.show()
