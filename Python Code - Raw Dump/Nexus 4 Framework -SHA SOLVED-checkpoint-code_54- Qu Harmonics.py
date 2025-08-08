import hashlib

def compute_sha512(input_str):
    """Compute the SHA-512 hash for the given input string."""
    return hashlib.sha512(input_str.encode()).hexdigest()

def similarity_score(hash1, hash2):
    """Calculate the percentage similarity between two hashes."""
    return sum(1 for a, b in zip(hash1, hash2) if a == b) / len(hash1) * 1000

def refine_input(seed_input, target_hash, max_iterations=100):
    """Refine the input string iteratively to match the target hash."""
    current_input = seed_input
    best_similarity = 0
    best_input = current_input

    for iteration in range(1, max_iterations + 1):
        # Generate a mutated candidate input
        candidate_input = mutate_input(current_input)
        
        # Compute hash and similarity score
        candidate_hash = compute_sha512(candidate_input)
        similarity = similarity_score(candidate_hash, target_hash)

        print(f"Iteration {iteration}: Similarity = {similarity:.2f}%")
        
        if similarity > best_similarity:
            best_similarity = similarity
            best_input = candidate_input

            # Stop if a perfect match is found
            if similarity == 100.0:
                print(f"Converged after {iteration} iterations!")
                break

        # Update the input for the next iteration
        current_input = candidate_input

    return best_input, best_similarity

def mutate_input(input_str):
    """Mutate the input string by randomly altering one character."""
    import random
    mutation_index = random.randint(0, len(input_str) - 1)
    mutated_char = chr(random.randint(32, 126))  # Replace with a random printable character
    return input_str[:mutation_index] + mutated_char + input_str[mutation_index + 1:]

# Example usage
seed_input = "Watson, the game is afoot"
target_hash = "ee1f36c84231330aa9ea6a842a36ff2eb09d5d004c0f449e11df865f309a3ac9d71f115b981eda6d4cfaa531434e4e5f4a73a61720430e8763700e14a1387a05"

final_input, final_similarity = refine_input(seed_input, target_hash)
print(f"Final Approximated Input: {final_input}")
print(f"Final Similarity: {final_similarity:.2f}%")
