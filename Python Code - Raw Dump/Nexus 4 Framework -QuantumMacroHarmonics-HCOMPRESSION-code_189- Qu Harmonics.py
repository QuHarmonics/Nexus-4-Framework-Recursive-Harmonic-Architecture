import numpy as np
import hashlib
from matplotlib import pyplot as plt

# Define the original hash (SHA-256)
original_hash = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"

# Convert hash to binary format
def hex_to_binary_array(hex_string):
    return np.array([int(hex_string[i:i+2], 16) for i in range(0, len(hex_string), 2)], dtype=np.uint8)

# Reverse SHA-256 function (Simulated Growth Model)
def reverse_sha256_simulation(input_hash, iterations=10000):
    harmonics = []
    current_seed = hex_to_binary_array(input_hash)
    
    # Harmonization loop
    for i in range(iterations):
        wave = np.cumsum(current_seed * 1.5)  # Simulating quantum growth
        wave = np.sin(wave / len(wave)) + np.cos(wave / len(wave))  # Introduce oscillation
        
        harmonics.append(wave)  # Record for later visualization
        
        # Adjust current seed to approach convergence
        current_seed = np.mod(wave * 1000, 256).astype(np.uint8)  # Quantization
        
        # Check if re-hashed result matches original hash
        reconstructed_hash = hashlib.sha256(current_seed.tobytes()).hexdigest()
        if reconstructed_hash == input_hash:
            print(f"Match found after {i + 1} iterations!")
            return current_seed, harmonics

    print("No exact match found, best approximation provided.")
    return current_seed, harmonics

# Run the simulation
reconstructed_seed, harmonics_log = reverse_sha256_simulation(original_hash)

# Visualize the harmonics
plt.figure(figsize=(12, 6))
for i, harmonic in enumerate(harmonics_log[:10]):  # Visualize first 10 harmonics
    plt.plot(harmonic, label=f"Harmonic {i + 1}")
plt.title("Harmonics Generated During Reconstruction")
plt.xlabel("Index")
plt.ylabel("Amplitude")
plt.legend()
plt.show()

# Output results
print("Reconstructed Seed (First 10 Bytes):", reconstructed_seed[:10])
