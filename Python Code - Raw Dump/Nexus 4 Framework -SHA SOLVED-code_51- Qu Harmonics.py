import hashlib
import random

def generate_hash(input_value):
    """Generate a SHA-512 hash."""
    return hashlib.sha512(input_value.encode('utf-8')).hexdigest()

def harmonic_refine_input(target_hash, initial_input, max_iterations=100, tolerance=0.001):
    """
    Refine input using harmonic constraints and feedback to approximate the target hash.
    """
    def similarity_score(hash1, hash2):
        # Compare hashes and calculate similarity
        return sum(c1 == c2 for c1, c2 in zip(hash1, hash2)) / len(hash1)

    current_input = initial_input
    iteration = 0
    history = []

    while iteration < max_iterations:
        # Hash the current input
        current_hash = generate_hash(current_input)
        history.append((current_input, current_hash))
        
        # Calculate similarity to the target hash
        similarity = similarity_score(current_hash, target_hash)
        print(f"Iteration {iteration + 1}: Similarity = {similarity:.2%}")

        # Check for convergence
        if similarity >= 1 - tolerance:
            print(f"Converged after {iteration + 1} iterations!")
            return current_input, history

        # Refine input using harmonic feedback
        # Modify input harmonically (e.g., alter characters based on structured patterns)
        input_list = list(current_input)
        for i in range(len(input_list)):
            # Adjust character harmonically (e.g., cycle through ASCII range)
            input_list[i] = chr((ord(input_list[i]) + random.randint(1, 3)) % 128)
        
        current_input = ''.join(input_list)
        iteration += 1

    return current_input, history

# Example usage
target_hash = "ee1f36c84231330aa9ea6a842a36ff2eb09d5d004c0f449e11df865f309a3ac9d71f115b981eda6d4cfaa531434e4e5f4a73a61720430e8763700e14a1387a05"
initial_input = "3333333333333333"  # Structured initial guess

# Perform harmonic refinement
approx_input, process_history = harmonic_refine_input(target_hash, initial_input)

print("Final Approximated Input:", approx_input)
