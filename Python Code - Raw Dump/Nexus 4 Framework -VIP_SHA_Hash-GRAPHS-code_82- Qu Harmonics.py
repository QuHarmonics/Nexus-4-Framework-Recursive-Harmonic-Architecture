import math

def recursive_update(H_prev, n, T, alpha, gamma):
    # Calculate the error and the gain factor based on that error.
    error = T - H_prev
    gain = 1 + gamma * abs(error)
    
    # The oscillatory correction: note that cos(n*pi) alternates between 1 and -1.
    oscillatory_term = H_prev * (-0.5) * math.cos(n * math.pi)
    
    # The error correction term, damped by 1/(n+1) and scaled by gain.
    correction_term = alpha * (error / (n + 1)) * gain
    
    # The new state is the sum of the oscillatory term and the correction term.
    return oscillatory_term + correction_term

# Parameters
T = 0.35         # Target harmonic constant
alpha = 1.5      # Amplification factor
gamma = 1.0      # Gain scaling factor
iterations = 5000  # Total number of iterations to simulate
H = 0.5          # Initial state (starting above the target to test convergence)

# List to store the sequence of H values
results = [H]

# Run the recursive update for the specified number of iterations
for n in range(1, iterations + 1):
    H = recursive_update(H, n, T, alpha, gamma)
    results.append(H)
    print(f"Iteration {n}: H = {H:.8f}")

print("\nFinal state after", iterations, "iterations:", H)
