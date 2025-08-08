import numpy as np
import hashlib

# Hash Input Functions
def sha512_hash(data):
    return np.frombuffer(hashlib.sha512(data.encode('utf-8')).digest(), dtype=np.uint8)

def large_data_hash(data):
    # Use SHA-512 on chunks of large data to create a combined hash
    chunk_size = 1024  # Example chunk size
    combined_hash = hashlib.sha512()
    for i in range(0, len(data), chunk_size):
        combined_hash.update(data[i:i + chunk_size].encode('utf-8'))
    return np.frombuffer(combined_hash.digest(), dtype=np.uint8)

# Samson-based Outer Loop
def samson_outer_loop(small_hash, large_hash, max_iterations=100):
    best_params = None
    best_score = float('inf')  # Lower score is better
    meta_params = np.linspace(0.5, 2.0, 10)  # Example parameter range

    for param in meta_params:
        print(f"Testing meta-parameter: {param}")

        # Run the inner loop with the current parameter
        small_result, small_score = inner_loop(small_hash, param)
        large_result, large_score = inner_loop(large_hash, param)

        # Evaluate combined score
        combined_score = small_score + large_score
        print(f"Combined score: {combined_score}")

        # Check for improvement
        if combined_score < best_score:
            best_score = combined_score
            best_params = param

    print(f"Optimal meta-parameter found: {best_params}")
    return best_params

# Inner Loop (Unfolding Logic)
def inner_loop(hash_input, meta_param):
    # Ensure hash_input is a numeric array
    hash_input = hash_input.astype(np.float64)

    harmonics = np.cumsum(hash_input * meta_param)
    reconstructed = np.round(np.diff(harmonics) / meta_param).astype(np.uint8)
    reconstructed = np.insert(reconstructed, 0, harmonics[0] / meta_param)

    # Calculate error or stability score
    error = np.sum(np.abs(hash_input - reconstructed))
    return reconstructed, error

# Generate input hashes
small_data = "This is a small dataset."
large_data = "This is a very large dataset " * 1000  # Large dataset
small_hash = sha512_hash(small_data)
large_hash = large_data_hash(large_data)

# Run Samson-guided tuning
optimal_param = samson_outer_loop(small_hash, large_hash)

# Final run with the optimal parameter
final_small_result, _ = inner_loop(small_hash, optimal_param)
final_large_result, _ = inner_loop(large_hash, optimal_param)

# Display results
print("Final Reconstruction (Small Hash):", final_small_result)
print("Final Reconstruction (Large Hash):", final_large_result)
