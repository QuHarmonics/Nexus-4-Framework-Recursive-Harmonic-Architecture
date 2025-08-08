import hashlib

def text_to_hex(text):
    """
    Converts text to hexadecimal.
    """
    return ''.join(format(ord(c), '02x') for c in text)

def hex_to_text(hex_str):
    """
    Converts hexadecimal to text.
    """
    try:
        return bytes.fromhex(hex_str).decode('utf-8')
    except ValueError:
        return hex_str  # Keep as-is if conversion fails

def recursive_tuning(input_value, guide_value, max_iterations=100, tolerance=0.01):
    """
    Recursively tunes the input value to align with the guide value.
    - input_value: Starting value for the tuning machine.
    - guide_value: The pre-hash guide value (target).
    - max_iterations: Maximum number of iterations to run the loop.
    - tolerance: Threshold for stopping based on similarity.
    """
    def similarity_score(a, b):
        """
        Calculates a similarity score between two strings as a percentage.
        """
        min_len = min(len(a), len(b))
        score = sum(1 for i in range(min_len) if a[i] == b[i]) / min_len
        return score

    current_value = input_value
    iteration = 0
    history = [current_value]  # To track progress

    while iteration < max_iterations:
        # Convert to hex and back, simulating harmonic unfolding
        current_value = text_to_hex(current_value)
        current_value = hex_to_text(current_value)

        # Compare with the guide value
        score = similarity_score(current_value, guide_value)

        print(f"Iteration {iteration + 1}: Similarity = {score:.2%}")
        print(f"Current Value: {current_value}")

        history.append(current_value)

        # Check if similarity meets tolerance
        if score >= 1.0 - tolerance:
            print(f"Alignment achieved at iteration {iteration + 1}!")
            break

        iteration += 1

    return current_value, history

# Example inputs
pre_hash_512 = "33333332333333353334333633333335333333363333333133343331"  # Guide value (pre-hash)
input_value = "33333332333333353334333633333335333333363333333133343331"  # Starting value

# Run the tuning machine
final_result, process_history = recursive_tuning(input_value, pre_hash_512, max_iterations=100, tolerance=0.001)

# Output results
print("\nFinal Harmonized Value:", final_result)
print("\nProcess History:")
for i, step in enumerate(process_history):
    print(f"Step {i}: {step}")
