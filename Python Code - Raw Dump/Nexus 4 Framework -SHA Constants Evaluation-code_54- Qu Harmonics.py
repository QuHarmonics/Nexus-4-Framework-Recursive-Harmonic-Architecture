import hashlib
import numpy as np
import matplotlib.pyplot as plt

# Generate a real SHA-256 hash from a known input
input_data = "hello"
real_hash = hashlib.sha256(input_data.encode()).hexdigest()
print(f"Real Hash (from '{input_data}'): {real_hash}")

# Constants and setup
harmonic_constant = 0.35
iterations = 10  # Integral triangle bound
reverse_waveform = np.random.rand(8, 8)  # Randomized starting lattice
constants = np.array([0.27264203, 0.46389402, 0.74472339, 0.9576116,
                      0.23494206, 0.36852961, 0.59924109, 0.7011437])

# Recursive feedback refinement with real hash
entropy_history = []
for i in range(iterations):
    # Harmonic adjustment
    reverse_waveform += harmonic_constant * np.sin(reverse_waveform + constants[i % len(constants)])
    reverse_waveform = np.mod(reverse_waveform, 1.0)  # Ensure values stay in range

    # Calculate entropy (std deviation)
    entropy = np.std(reverse_waveform)
    entropy_history.append(entropy)

    # Compute the hash of the current waveform
    current_hash = hashlib.sha256(reverse_waveform.tobytes()).hexdigest()

    # Stop if hash matches
    if current_hash == real_hash:
        print(f"Hash matched at iteration {i + 1}")
        break

# Visualize entropy trends
plt.plot(range(1, len(entropy_history) + 1), entropy_history, marker='o')
plt.axhline(harmonic_constant, color='red', linestyle='--', label="Harmonic Target (0.35)")
plt.title("Entropy Convergence with Real Hash")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid(True)
plt.show()

# Output results
print("Final Reconstructed Hash:", current_hash)
print("Real Hash:", real_hash)
print("Final Waveform:")
print(reverse_waveform)
