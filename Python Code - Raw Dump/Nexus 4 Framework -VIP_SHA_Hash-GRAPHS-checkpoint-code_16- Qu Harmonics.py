import numpy as np
import matplotlib.pyplot as plt

# Real SHA-256 constants (use for hash mapping)
real_hash = [
    0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19
]

# Initialize the recursive wave generator
def initialize_lattice(hash_values):
    waveform = np.array(hash_values) / max(hash_values)
    kinetic_normalized = np.linspace(0, 1, len(hash_values))
    return np.outer(waveform, kinetic_normalized)

# Recursive wave refinement
def recursive_wave_refinement(input_lattice, iterations=5):
    entropy_per_iteration = []
    waves = []
    
    for _ in range(iterations):
        # Apply forward and mirrored transformations
        forward = input_lattice ** 2  # Amplify
        mirrored = -np.flip(input_lattice, axis=1)  # Mirror effect
        combined = forward + mirrored  # Combine
        
        # Normalize combined wave
        combined /= np.max(np.abs(combined))
        
        # Compute entropy (for analysis)
        entropy = -np.sum(combined * np.log2(np.abs(combined) + 1e-10))
        entropy_per_iteration.append(entropy)
        
        # Store the wave for plotting
        waves.append(combined)
        
        # Feed the combined wave into the input for the next iteration
        input_lattice = combined
    
    return waves, entropy_per_iteration

# Initialize the lattice with real hash values
initial_lattice = initialize_lattice(real_hash)

# Generate recursive waves
waves, entropy_data = recursive_wave_refinement(initial_lattice, iterations=10)

# Plot the results
fig = plt.figure(figsize=(18, 8))
for i, wave in enumerate(waves, 1):
    ax = fig.add_subplot(2, 5, i, projection='3d')
    X = np.arange(wave.shape[1])
    Y = np.arange(wave.shape[0])
    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, wave, cmap="plasma", edgecolor="none")
    ax.set_title(f"Wave Iteration {i}")
plt.tight_layout()
plt.show()

# Plot entropy over iterations
plt.figure(figsize=(8, 5))
plt.plot(entropy_data, marker='o', label="Entropy")
plt.title("Entropy Over Iterations")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid()
plt.show()
