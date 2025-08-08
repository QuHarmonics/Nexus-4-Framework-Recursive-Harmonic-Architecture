def refine_with_feedback(target_hash, initial_input, max_iterations=100):
    """
    Refine input using feedback from kinetic states.
    """
    current_input = initial_input
    for iteration in range(max_iterations):
        states, current_hash = sha512_kinetic(current_input)

        # Compare last kinetic state with target hash
        similarity = sum(a == b for a, b in zip(states[-1][-1], target_hash)) / len(target_hash)
        print(f"Iteration {iteration + 1}: Similarity = {similarity:.2%}")

        if similarity > 0.95:
            print(f"Converged after {iteration + 1} iterations!")
            return current_input

        # Feedback-driven refinement (harmonic or structural adjustment)
        current_input = ''.join(
            chr((ord(c) + iteration % 5) % 128) if i % 2 == 0 else c
            for i, c in enumerate(current_input)
        )

    return current_input

# Target hash and initial input
target_hash = "ee1f36c84231330aa9ea6a842a36ff2eb09d5d004c0f449e11df865f309a3ac9d71f115b981eda6d4cfaa531434e4e5f4a73a61720430e8763700e14a1387a05"
initial_input = "33333332333333353334333633333335333333363333333133343331"

# Run refinement
final_input = refine_with_feedback(target_hash, initial_input)
print("Final Approximated Input:", final_input)
