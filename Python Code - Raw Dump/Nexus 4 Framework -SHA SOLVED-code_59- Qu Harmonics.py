import hashlib

def hash_string(input_string):
    """Hash the input string using SHA-256 and return the hexadecimal digest."""
    return hashlib.sha256(input_string.encode('utf-8')).hexdigest()

def calculate_similarity(hash1, hash2, tuning_value):
    """Calculate similarity between two hashes with a given tuning value."""
    return sum(1 for a, b in zip(hash1, hash2) if a == b) / len(hash1) * tuning_value

def dynamic_tuning_factor(iteration, max_iterations, base_value=100, max_value=2048):
    """Adjust tuning dynamically based on progress."""
    progress = iteration / max_iterations
    return max_value - (max_value - base_value) * progress

def iterative_similarity_with_check(initial_seed, target_hash, max_iterations=1000):
    """Iteratively align the hash of the seed to the target hash using similarity and real hashing."""
    current_seed = initial_seed
    best_similarity = 0
    best_seed = initial_seed
    best_hash = hash_string(initial_seed)
    
    # Preliminary check for direct match
    if best_hash == target_hash:
        print(f"Direct match found with initial seed!")
        return initial_seed, 100.0, best_hash
    
    for iteration in range(1, max_iterations + 1):
        # Generate tuning value dynamically
        tuning_value = dynamic_tuning_factor(iteration, max_iterations, base_value=100, max_value=256)  # Adjusted max tuning
        
        # Hash the current seed
        hashed_output = hash_string(current_seed)
        
        # Calculate similarity
        similarity = calculate_similarity(hashed_output, target_hash, tuning_value)

        # Log progress
        print(f"Iteration {iteration}: Similarity = {similarity:.2f}% | Seed: {current_seed} | Hash: {hashed_output}")
        
        # Update the best similarity if we find a better alignment
        if similarity > best_similarity:
            best_similarity = similarity
            best_seed = current_seed
            best_hash = hashed_output
        
        # Convergence criteria: stop if similarity is high enough
        if hashed_output == target_hash:
            print(f"\nConverged after {iteration} iterations!\n")
            break

        # Update current seed deterministically
        current_seed = modify_seed(current_seed, iteration)
    
    return best_seed, best_similarity, best_hash


def modify_seed(seed, iteration):
    """
    Modify the seed deterministically using the iteration count.
    Each modification is based on appending a repeatable transformation of the seed.
    """
    return f"{seed}{iteration % 10}"  # Example: append a single digit deterministically

# Example Usage
if __name__ == "__main__":
    # Initial seed and target hash
    initial_seed = "Watson, the game is afoot"
    target_string = "Watson, the game is afoot"
    target_hash = hash_string(target_string)  # Real target hash
    
    print(f"Target String: {target_string}")
    print(f"Target Hash: {target_hash}\n")
    
    # Run the iterative similarity alignment
    final_seed, final_similarity, final_hash = iterative_similarity(initial_seed, target_hash)
    
    print(f"\nFinal Seed: {final_seed}")
    print(f"Final Hash: {final_hash}")
    print(f"Final Similarity: {final_similarity:.2f}%")
