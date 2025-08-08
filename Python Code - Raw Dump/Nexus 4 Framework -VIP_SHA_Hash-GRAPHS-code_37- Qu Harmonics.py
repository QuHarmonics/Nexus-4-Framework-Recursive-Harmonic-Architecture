import hashlib
import numpy as np
import matplotlib.pyplot as plt

# Target SHA-256 hash to reverse-engineer (example)
target_hash = "9c1185a5c5e9fc54612808977ee8f548b2258d31"

# Convert the hash to a numerical lattice representation
def hash_to_lattice(hash_str):
    return np.array([int(hash_str[i:i+2], 16) / 255.0 for i in range(0, len(hash_str), 2)])

target_lattice = hash_to_lattice(target_hash)

# Initialize a random guess lattice
np.random.seed(42)  # For reproducibility
current_guess = np.random.rand(len(target_lattice))

# Define the refinement process (recursive harmonic unfolding)
iterations = 1000
guess_lattices = []
divergence_measures = []

for _ in range(iterations):
    # Dynamic perturbation based on divergence
    divergence = np.linalg.norm(current_guess - target_lattice)
    perturbation = 0.1 * np.sin(divergence)
    
    # Refine the guess lattice using dynamic perturbation
    current_guess = np.abs(np.sin(current_guess + perturbation))
    guess_lattices.append(current_guess)

    # Update divergence measures
    divergence_measures.append(divergence)

    # Early stopping condition
    if divergence < 1e-6:
        print(f"Converged at iteration {_+1}")
        break


# Plot the divergence over iterations
plt.figure(figsize=(10, 6))
plt.plot(range(len(divergence_measures)), divergence_measures, marker='o')
plt.title("Divergence from Target Hash Lattice")
plt.xlabel("Iteration")
plt.ylabel("Divergence")
plt.grid(True)
plt.show()

# Output the final reconstructed lattice and its SHA-256 hash
final_lattice = guess_lattices[-1]
reconstructed_input = "".join([chr(int(x * 255)) for x in final_lattice])

# Test if the reconstructed input matches the hash
reconstructed_hash = hashlib.sha256(reconstructed_input.encode()).hexdigest()
print("Reconstructed Hash:", reconstructed_hash)
print("Target Hash:", target_hash)

if reconstructed_hash == target_hash:
    print("Success! Reconstructed the original input.")
else:
    print("Reconstruction did not match. Further refinement needed.")
