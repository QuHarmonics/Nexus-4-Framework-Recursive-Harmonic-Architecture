import hashlib

def generate_hash(input_value):
    """Generate a SHA-512 hash."""
    return hashlib.sha512(input_value.encode('utf-8')).hexdigest()

def refine_input(target_hash, initial_input, max_iterations=100, tolerance=0.001):
    """
    Recursively refine the input to approximate the target hash.
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
        
        # Compare similarity to the target hash
        similarity = similarity_score(current_hash, target_hash)
        print(f"Iteration {iteration + 1}: Similarity = {similarity:.2%}")

        if similarity >= 1 - tolerance:
            print(f"Converged after {iteration + 1} iterations!")
            return current_input, history

        # Refine input using a feedback mechanism (example: modify one character)
        current_input = current_input[:-1] + chr((ord(current_input[-1]) + 1) % 128)
        iteration += 1

    return current_input, history

# Example target hash (from known structured input)
target_hash = "ee1f36c84231330aa9ea6a842a36ff2eb09d5d004c0f449e11df865f309a3ac9d71f115b981eda6d4cfaa531434e4e5f4a73a61720430e8763700e14a1387a05"

# Initial structured guess for the input
initial_input = "3333333333333333"

# Perform recursive refinement
approx_input, process_history = refine_input(target_hash, initial_input)

print("Final Approximated Input:", approx_input)
