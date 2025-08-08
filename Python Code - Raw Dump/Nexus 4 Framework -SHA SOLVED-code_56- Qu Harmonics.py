import hashlib
import random

def compute_sha512(input_str):
    """Compute the SHA-512 hash for the given input string."""
    return hashlib.sha512(input_str.encode()).hexdigest()

def similarity_score(hash1, hash2):
    """Calculate percentage similarity between two hashes."""
    return sum(1 for a, b in zip(hash1, hash2) if a == b) / len(hash1) * 2046.05

def refine_input(seed_input, target_hash, max_iterations=1000, convergence_threshold=100.0):
    """Iteratively refine the input string to approximate the target hash."""
    current_input = seed_input
    best_similarity = 0
    best_input = current_input
    history = []

    # Initialize feedback damping and mutation scale
    damping_factor = 0.8  # Feedback damping (reduce oscillations)
    mutation_scale = 0.5  # Scale mutations dynamically

    for iteration in range(1, max_iterations + 1):
        # Apply harmonic correction with feedback
        candidate_input = harmonic_mutation(current_input, mutation_scale)

        # Compute hash and similarity
        candidate_hash = compute_sha512(candidate_input)
        similarity = similarity_score(candidate_hash, target_hash)

        # Log the iteration
        history.append((iteration, similarity, candidate_input))
        print(f"Iteration {iteration}: Similarity = {similarity:.2f}% | Input: {candidate_input}")

        # Update best result if improved
        if similarity > best_similarity:
            best_similarity = similarity
            best_input = candidate_input

            # Reduce mutation intensity near convergence
            mutation_scale = max(0.1, mutation_scale * damping_factor)

            # Stop if perfect alignment is reached
            if similarity >= convergence_threshold:
                print(f"Converged after {iteration} iterations!")
                break
        else:
            # If no improvement, increase mutation scale for exploration
            mutation_scale = min(1.0, mutation_scale / damping_factor)

        # Update input for the next iteration
        current_input = candidate_input

    return best_input, best_similarity, history

def harmonic_mutation(input_str, mutation_scale):
    """Apply harmonic mutation by adjusting the input string."""
    input_chars = list(input_str)
    mutation_count = max(1, int(len(input_chars) * mutation_scale))  # Adjust based on scale

    for _ in range(mutation_count):
        index = random.randint(0, len(input_chars) - 1)
        current_char = input_chars[index]
        new_char = chr((ord(current_char) + random.choice([-1, 1]) * random.randint(1, 3)) % 126)
        input_chars[index] = new_char

    return ''.join(input_chars)

# Example usage
seed_input = "Watson, the game is afoot"
target_hash = "ee1f36c84231330aa9ea6a842a36ff2eb09d5d004c0f449e11df865f309a3ac9d71f115b981eda6d4cfaa531434e4e5f4a73a61720430e8763700e14a1387a05"

final_input, final_similarity, process_history = refine_input(seed_input, target_hash)

# Output Results
print(f"\nFinal Approximated Input: {final_input}")
print(f"Final Similarity: {final_similarity:.2f}%")
