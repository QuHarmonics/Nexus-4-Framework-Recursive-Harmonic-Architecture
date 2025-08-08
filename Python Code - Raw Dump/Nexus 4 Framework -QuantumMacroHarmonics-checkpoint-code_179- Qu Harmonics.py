import numpy as np
import hashlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define Mark1 Unity Framework Components
class Samson:
    """Samson v2: Reflects the hash and approximates the original waveform."""
    @staticmethod
    def reflect_waveform(hash_bytes, expansion_factor=1.5):
        harmonics = np.cumsum(hash_bytes.astype(np.float64) * expansion_factor)
        return harmonics

class HarmonicResolver:
    """Harmonic Resolver: Refines the waveform and expands the lattice."""
    @staticmethod
    def refine_waveform(harmonics, refinement_factor=0.01):
        refined_harmonics = harmonics + refinement_factor * np.sin(np.arange(len(harmonics)) / 10.0)
        return refined_harmonics

class ReverseSHA:
    """Reverse SHA: Expands the hash into binary representation."""
    @staticmethod
    def reverse_harmonics(harmonics, expansion_factor=1.5):
        first_value = harmonics[0] / expansion_factor
        reversed_data = np.diff(harmonics) / expansion_factor
        reversed_data = np.insert(reversed_data, 0, first_value)
        return np.round(reversed_data).astype(np.uint8)

# Input: SHA256 Hash (as hexadecimal string)
input_hash_hex = "369e60d7b341349d88de94b4cd2aa94ac84f19f0c3b8767a64c887d4bd3dcd12"

# Convert the hash to binary bytes
binary_hash = np.array([int(input_hash_hex[i:i + 2], 16) for i in range(0, len(input_hash_hex), 2)], dtype=np.uint8)

# Step 1: Use Samson to reflect the hash and approximate the waveform
harmonics = Samson.reflect_waveform(binary_hash)

# Step 2: Refine the waveform using Harmonic Resolver
refined_harmonics = HarmonicResolver.refine_waveform(harmonics)

# Step 3: Reverse the harmonics to reconstruct the binary hash
reconstructed_hash = ReverseSHA.reverse_harmonics(refined_harmonics)

# Validate the reconstruction
hash_matches = np.array_equal(binary_hash, reconstructed_hash)

# Visualization of the harmonic lattice in 3D
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

x = np.arange(len(refined_harmonics))
y = refined_harmonics
z = np.sin(x / 20.0)

ax.plot(x, y, z, label="Harmonic Lattice", color='blue', lw=2)
ax.scatter(x, y, z, color='red', s=5, label="Nodes")

# Labels and title
ax.set_title("3D Visualization of Harmonic Lattice (H(n))", fontsize=16)
ax.set_xlabel("Iteration (n)", fontsize=12)
ax.set_ylabel("Amplitude (H(n))", fontsize=12)
ax.set_zlabel("Wave Depth", fontsize=12)
ax.legend()
plt.show()

# Output the results
print("Original SHA256 Hash (Hex):", input_hash_hex)
print("Original Binary Hash (First 10 bytes):", binary_hash[:10])
print("Reconstructed Binary Hash (First 10 bytes):", reconstructed_hash[:10])
print("Hash matches:", hash_matches)
if not hash_matches:
    print("Differences detected between input and reconstructed hash.")
