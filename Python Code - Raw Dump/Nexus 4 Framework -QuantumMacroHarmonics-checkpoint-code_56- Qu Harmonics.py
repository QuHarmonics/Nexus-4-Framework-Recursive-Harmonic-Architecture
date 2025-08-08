import numpy as np
import matplotlib.pyplot as plt
from mpmath import zetazero

# Recursive realignment using Samson's reflective convergence
def samson_realign_waveform(waveform, target_ratio=0.35, tolerance=1e-6, max_iterations=1e6):
    """
    Realigns a waveform to achieve harmonic stability using Samson's recursive reflection.
    """
    iteration = 0
    previous_waveform = waveform.copy()
    current_waveform = waveform.copy()
    stability_ratio = float('inf')

    while stability_ratio > tolerance and iteration < max_iterations:
        # Dynamically calculate the zeta zeros for this iteration
        num_zeros = iteration + 1
        zeta_zeros = [zetazero(n).real for n in range(1, num_zeros + 1)]

        # Reflect the waveform using zeta zeros
        reflected_waveform = reflect_waveform_samson(current_waveform, zeta_zeros)

        # Calculate stability ratio (ΣPᵢ/ΣAᵢ ≈ 0.35)
        P = np.sum(np.abs(reflected_waveform - previous_waveform))
        A = np.sum(np.abs(previous_waveform))
        stability_ratio = P / A if A != 0 else float('inf')

        # Check if stability ratio is within the target range
        if abs(stability_ratio - target_ratio) < tolerance:
            break

        # Update waveforms for the next iteration
        previous_waveform = current_waveform.copy()
        current_waveform = reflected_waveform
        iteration += 1

    return current_waveform, stability_ratio, iteration

# Reflect waveform based on zeta zeros
def reflect_waveform_samson(waveform, zeta_zeros):
    reflected_waveform = waveform.copy()
    harmonic_length = len(waveform)
    for zero in zeta_zeros:
        zero_index = int((zero / max(zeta_zeros)) * harmonic_length)
        if 0 <= zero_index < harmonic_length:
            # Shift the waveform to reflect at the zero
            reflected_waveform[zero_index:] = np.roll(reflected_waveform[zero_index:], shift=1)
    return reflected_waveform

# Convert a hash into binary and then into a waveform
def hash_to_binary(hash_hex):
    return np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])

def binary_to_waveform(binary_data):
    n = len(binary_data)
    x = np.linspace(0, 2 * np.pi, n)
    return np.sin(x) * (binary_data - 0.5)

# Example: Realign the hash waveform
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"  # SHA-256 for 'abc'

# Generate the original waveform from hash
binary_data = hash_to_binary(hash_hex)
hash_waveform = binary_to_waveform(binary_data)

# Realign the waveform using Samson
realigned_waveform, final_stability_ratio, total_iterations = samson_realign_waveform(hash_waveform)

# Visualization of the realignment process
plt.figure(figsize=(12, 10))

plt.subplot(3, 1, 1)
plt.plot(hash_waveform[:128], label="Original Hash Waveform", color='orange')
plt.title("Original Hash Waveform (First 128 Points)")
plt.grid()
plt.legend()

plt.subplot(3, 1, 2)
plt.plot(realigned_waveform[:128], label="Realigned Waveform (Samson)", color='green')
plt.title(f"Realigned Waveform (First 128 Points) - Stability Ratio: {final_stability_ratio:.6f}")
plt.grid()
plt.legend()

plt.subplot(3, 1, 3)
plt.plot(np.abs(realigned_waveform[:128] - hash_waveform[:128]), label="Difference (Realigned vs. Original)", color='blue')
plt.title(f"Difference Waveform (First 128 Points) - Total Iterations: {total_iterations}")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()

# Output results
print("Final Stability Ratio:", final_stability_ratio)
print("Total Iterations:", total_iterations)
