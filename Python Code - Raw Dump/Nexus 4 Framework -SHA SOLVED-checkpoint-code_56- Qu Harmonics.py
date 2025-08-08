import random

def dynamic_tuning_factor(iteration, max_iterations, base_value=100, max_value=2048):
    """Adjust tuning dynamically based on progress."""
    progress = iteration / max_iterations
    # Gradually reduce from max_value to base_value
    return max_value - (max_value - base_value) * progress

def calculate_similarity(hash1, hash2, tuning_value):
    """Calculate similarity between two hashes with a given tuning value."""
    return sum(1 for a, b in zip(hash1, hash2) if a == b) / len(hash1) * tuning_value

def generate_random_hash(length=32):
    """Generate a random hash-like string."""
    return ''.join(random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") for _ in range(length))

def iterative_similarity(hash1, target_hash, max_iterations=1000):
    """Iteratively align hash1 to target_hash using similarity tuning."""
    current_hash = hash1
    best_similarity = 0
    best_hash = hash1
    
    for iteration in range(1, max_iterations + 1):
        # Generate tuning value dynamically
        tuning_value = dynamic_tuning_factor(iteration, max_iterations)

        # Mutate current hash slightly to explore new states
        mutated_hash = mutate_hash(current_hash)

        # Calculate similarity
        similarity = calculate_similarity(mutated_hash, target_hash, tuning_value)

        # Log progress
        print(f"Iteration {iteration}: Similarity = {similarity:.2f}% | Input: {mutated_hash}")

        # Update the best similarity if we find a better alignment
        if similarity > best_similarity:
            best_similarity = similarity
            best_hash = mutated_hash
        
        # Convergence criteria: stop if similarity is 100%
        if similarity >= 100.0:
            print(f"\nConverged after {iteration} iterations!\n")
            break

        # Update current hash for next iteration
        current_hash = mutated_hash

    return best_hash, best_similarity

def mutate_hash(hash_string):
    """Randomly mutate a single character in the hash string."""
    hash_list = list(hash_string)
    index = random.randint(0, len(hash_list) - 1)
    hash_list[index] = random.choice("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789")
    return ''.join(hash_list)

# Example Usage
if __name__ == "__main__":
    # Initialize hashes
    initial_hash = generate_random_hash()
    target_hash = "Watson,thegameisafoot"  # Target to align to
    
    print(f"Initial Hash: {initial_hash}")
    print(f"Target Hash: {target_hash}\n")
    
    # Run the iterative similarity alignment
    final_hash, final_similarity = iterative_similarity(initial_hash, target_hash)

    print(f"\nFinal Approximated Input: {final_hash}")
    print(f"Final Similarity: {final_similarity:.2f}%")
