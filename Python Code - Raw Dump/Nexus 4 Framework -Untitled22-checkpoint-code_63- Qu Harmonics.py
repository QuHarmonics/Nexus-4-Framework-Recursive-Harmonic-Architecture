# Correcting the dimensional mismatch and refining the lattice processing

def backpropagate_amplitude_phase_corrected(hash_val, lattice, iterations=10):
    # Converting hash value into a numerical representation
    hash_numeric = np.array([int(hash_val[i:i + 2], 16) for i in range(0, len(hash_val), 2)])
    lattice_size = lattice.shape[0]  # Assuming square lattice
    
    # Reshape to match lattice size by padding/truncating appropriately
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

# Run corrected unhashing process
def unhash_sha256_corrected(target_hash, size, iterations=10):
    lattice = initialize_reverse_lattice(size)
    lattice, amplitude_diff_history, phase_diff_history = backpropagate_amplitude_phase_corrected(target_hash, lattice, iterations)
    original_input = reconstruct_input(target_hash, lattice)
    return lattice, original_input, amplitude_diff_history, phase_diff_history

# Attempt with corrected implementation
size = 256  # For simplicity
iterations = 256

# Re-run the unhashing process
lattice, reconstructed_input, amplitude_diff_history, phase_diff_history = unhash_sha256_corrected(target_hash, size, iterations)

# Validate
is_valid = validate_reconstruction(reconstructed_input, target_hash)

# Plotting results
fig, axes = plt.subplots(3, 5, figsize=(20, 10))
for i in range(5):
    axes[0, i].imshow(amplitude_diff_history[i], cmap='viridis')
    axes[0, i].set_title(f'Amplitude Diff Iter {i+1}')
    axes[1, i].imshow(phase_diff_history[i], cmap='plasma')
    axes[1, i].set_title(f'Phase Diff Iter {i+1}')
    axes[2, i].imshow(lattice, cmap='inferno')
    axes[2, i].set_title(f'Lattice State Iter {i+1}')

plt.tight_layout()
plt.show()

# Display reconstructed input and validation result
reconstructed_input, is_valid
