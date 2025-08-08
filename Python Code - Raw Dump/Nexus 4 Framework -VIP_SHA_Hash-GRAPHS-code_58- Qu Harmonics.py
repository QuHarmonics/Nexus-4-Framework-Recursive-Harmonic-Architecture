import numpy as np
import matplotlib.pyplot as plt
import random

# --- Utilities ---

def initialize_reverse_lattice(size):
    return np.random.rand(size, size) * 255  # Values 0-255

def reconstruct_input(target_hash, lattice):
    flat_lattice = lattice.flatten()
    reconstructed_bytes = np.clip(np.round(flat_lattice), 0, 255).astype(np.uint8)
    reconstructed_hex = ''.join([format(b, '02x') for b in reconstructed_bytes[:len(target_hash)//2]])
    return reconstructed_hex

def validate_reconstruction(reconstructed_input, target_hash):
    return reconstructed_input.lower() == target_hash.lower()

# --- Core Correction ---

def backpropagate_amplitude_phase_corrected(hash_val, lattice, iterations=10):
    hash_numeric = np.array([int(hash_val[i:i + 2], 16) for i in range(0, len(hash_val), 2)])
    lattice_size = lattice.shape[0]
    
    if len(hash_numeric) < lattice_size**2:
        hash_numeric = np.pad(hash_numeric, (0, lattice_size**2 - len(hash_numeric)), mode='constant')
    elif len(hash_numeric) > lattice_size**2:
        hash_numeric = hash_numeric[:lattice_size**2]
    
    hash_numeric = hash_numeric.reshape((lattice_size, lattice_size))

    amplitude_diff_history = []
    phase_diff_history = []

    for _ in range(iterations):
        amplitude_diff = np.abs(lattice - hash_numeric)
        phase_diff = np.angle(np.exp(1j * (lattice - hash_numeric)))

        lattice -= 0.1 * amplitude_diff
        lattice += 0.1 * phase_diff

        amplitude_diff_history.append(np.copy(amplitude_diff))
        phase_diff_history.append(np.copy(phase_diff))

    return lattice, amplitude_diff_history, phase_diff_history

def unhash_sha256_corrected(target_hash, size, iterations=10):
    lattice = initialize_reverse_lattice(size)
    lattice, amplitude_diff_history, phase_diff_history = backpropagate_amplitude_phase_corrected(target_hash, lattice, iterations)
    original_input = reconstruct_input(target_hash, lattice)
    return lattice, original_input, amplitude_diff_history, phase_diff_history

# --- Simulate Target Hash ---

target_hash = ''.join(random.choices('0123456789abcdef', k=64))

# --- Run the Full Process ---

size = 16  # Keep small for visualization (256x256 would be huge for now)
iterations = 10

lattice, reconstructed_input, amplitude_diff_history, phase_diff_history = unhash_sha256_corrected(target_hash, size, iterations)

is_valid = validate_reconstruction(reconstructed_input, target_hash)

# --- Plotting Results ---

fig, axes = plt.subplots(3, 5, figsize=(20, 10))

for i in range(5):
    axes[0, i].imshow(amplitude_diff_history[i], cmap='viridis')
    axes[0, i].set_title(f'Amplitude Diff Iter {i+1}')
    axes[0, i].axis('off')
    
    axes[1, i].imshow(phase_diff_history[i], cmap='plasma')
    axes[1, i].set_title(f'Phase Diff Iter {i+1}')
    axes[1, i].axis('off')
    
    axes[2, i].imshow(lattice, cmap='inferno')
    axes[2, i].set_title(f'Lattice State After Iter {iterations}')
    axes[2, i].axis('off')

plt.tight_layout()
plt.show()

# --- Display Final Result ---

print("\nReconstructed Input:")
print(reconstructed_input)
print("\nValidation Result:")
print("Success!" if is_valid else "Failed")
