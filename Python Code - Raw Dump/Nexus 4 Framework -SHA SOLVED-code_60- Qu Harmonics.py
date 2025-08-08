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
