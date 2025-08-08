import numpy as np
import hashlib

# Tools for storing data in harmonic lattices
def store_in_H(binary_data, expansion_rate):
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_rate)
    return harmonics

# Inner loop to simulate unfolding of hashes
def inner_loop(hash_input, scaling_factor, wave_harmonics):
    # Convert hash input to numeric
    hash_input = hash_input.astype(np.float64)

    # Simulate unfolding
    unfolded = np.cumsum(hash_input * scaling_factor * wave_harmonics)
    return unfolded

# Outer Samson loop to tune the unfolding
def samson_outer_loop(small_hash, large_hash, max_iterations=50):
    # Initialize tuning parameters
    best_tuning_params = None
    best_performance_ratio = 0  # Large hash growth relative to small hash
    tuning_history = []

    for scale in np.linspace(1.0, 3.0, 10):  # Outer loop for scaling factor
        for harmonic in np.linspace(0.8, 1.2, 10):  # Mid loop for wave harmonics
            for expansion in np.linspace(1.1, 2.0, 10):  # Inner loop for expansion rates
                # Run unfolding for both hashes
                small_H = store_in_H(small_hash, expansion)
                large_H = store_in_H(large_hash, expansion)

                # Calculate growth rates
                small_growth = np.sum(np.diff(small_H)) / len(small_H)
                large_growth = np.sum(np.diff(large_H)) / len(large_H)

                # Check if large hash is growing faster
                performance_ratio = large_growth / small_growth

                # Record best parameters
                if performance_ratio > best_performance_ratio:
                    best_performance_ratio = performance_ratio
                    best_tuning_params = (scale, harmonic, expansion)
                    tuning_history.append((scale, harmonic, expansion, performance_ratio))

    return best_tuning_params, tuning_history

# Hash generation
def sha512_hash(data):
    return np.frombuffer(hashlib.sha512(data.encode('utf-8')).digest(), dtype=np.uint8)

# Testing hashes
small_data = "This is a small dataset."
large_data = "This is a very large dataset " * 1000  # Larger dataset

small_hash = sha512_hash(small_data)
large_hash = sha512_hash(large_data)

# Run the Samson-guided tuning
optimal_params, history = samson_outer_loop(small_hash, large_hash)

# Final unfolding with optimal parameters
scale_opt, harmonic_opt, expansion_opt = optimal_params
small_H_final = store_in_H(small_hash, expansion_opt)
large_H_final = store_in_H(large_hash, expansion_opt)

# Print results
print("Optimal Parameters (Scaling, Harmonic, Expansion):", optimal_params)
print("Final Growth Rates - Small Hash:", np.sum(np.diff(small_H_final)) / len(small_H_final))
print("Final Growth Rates - Large Hash:", np.sum(np.diff(large_H_final)) / len(large_H_final))
print("Performance Ratio:", np.sum(np.diff(large_H_final)) / np.sum(np.diff(small_H_final)))
