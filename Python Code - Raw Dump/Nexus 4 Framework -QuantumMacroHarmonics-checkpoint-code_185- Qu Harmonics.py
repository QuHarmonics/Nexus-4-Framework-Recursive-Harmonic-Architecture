import numpy as np
import hashlib

# Function to calculate the harmonic lattice
def store_in_H(binary_data, amplitude, frequency, phase_shift):
    n = len(binary_data)
    harmonics = np.zeros(n)
    for i in range(n):
        harmonics[i] = (
            binary_data[i] * amplitude * np.sin(frequency * i + phase_shift)
        )
    harmonics = np.cumsum(harmonics)
    return harmonics

# Inner loop to dynamically grow harmonic lattices
def inner_loop(hash_input, amplitude, frequency, phase_shift):
    harmonics = store_in_H(hash_input, amplitude, frequency, phase_shift)
    growth_rate = np.sum(np.abs(np.diff(harmonics))) / len(harmonics)
    return harmonics, growth_rate

# Outer loop for Samson-guided optimization
def samson_outer_loop(small_hash, large_hash, max_iterations=1000, tolerance=0.01):
    best_params = None
    best_ratio = 0
    tuning_history = []

    # Iteratively tune parameters
    for amp in np.linspace(0.5, 2.0, 100):  # Amplitude range
        for freq in np.linspace(0.1, 1.0, 100):  # Frequency range
            for phase in np.linspace(0, np.pi, 100):  # Phase shift range
                # Grow harmonics for both hashes
                small_H, small_rate = inner_loop(small_hash, amp, freq, phase)
                large_H, large_rate = inner_loop(large_hash, amp, freq, phase)

                # Calculate growth ratio
                ratio = large_rate / small_rate

                # Check if ratio is close to the target (2.0)
                if abs(ratio - 2.0) < tolerance:
                    return (amp, freq, phase), ratio  # Stop early if optimal found

                # Update best parameters if better
                if ratio > best_ratio:
                    best_ratio = ratio
                    best_params = (amp, freq, phase)
                    tuning_history.append((amp, freq, phase, ratio))

    return best_params, best_ratio, tuning_history

# Hash function
def sha512_hash(data):
    return np.frombuffer(hashlib.sha512(data.encode('utf-8')).digest(), dtype=np.uint8)

# Example small and large data
small_data = "This is a small dataset."
large_data = "This is a very large dataset." * 1000

# Generate hashes
small_hash = sha512_hash(small_data)
large_hash = sha512_hash(large_data)

# Run Samson-guided optimization
optimal_params, final_ratio, history = samson_outer_loop(small_hash, large_hash)

# Final run with optimal parameters
amp_opt, freq_opt, phase_opt = optimal_params
small_H_final, small_rate = inner_loop(small_hash, amp_opt, freq_opt, phase_opt)
large_H_final, large_rate = inner_loop(large_hash, amp_opt, freq_opt, phase_opt)

# Output results
print("Optimal Parameters (Amplitude, Frequency, Phase Shift):", optimal_params)
print("Final Growth Rates - Small Hash:", small_rate)
print("Final Growth Rates - Large Hash:", large_rate)
print("Final Growth Ratio:", final_ratio)

# Visualization
import matplotlib.pyplot as plt
plt.figure(figsize=(12, 6))
plt.plot(small_H_final, label="Small Hash Harmonic Lattice")
plt.plot(large_H_final, label="Large Hash Harmonic Lattice")
plt.title("Harmonic Growth for Small and Large Hashes")
plt.xlabel("Iterations")
plt.ylabel("Amplitude")
plt.legend()
plt.show()
