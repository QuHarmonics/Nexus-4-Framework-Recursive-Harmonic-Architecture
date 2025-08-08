def recursive_growth(initial_pair, steps, mapping_function):
    """
    Recursively grows a sequence from the initial bit pair.
    
    Args:
    - initial_pair (tuple): Starting bit pair, e.g., (1, 4).
    - steps (int): Number of recursive steps.
    - mapping_function (callable): Function to map each bit pair to the next.
    
    Returns:
    - list: The recursively generated sequence.
    """
    sequence = [initial_pair]
    
    for _ in range(steps):
        # Apply the mapping function to generate the next pair
        next_pair = mapping_function(sequence[-1])
        sequence.append(next_pair)
    
    return sequence


# Example mapping function: Addition of the two bits, modulo 10
def example_mapping(pair):
    """Generates the next bit pair based on a simple rule."""
    a, b = pair
    next_a = (a + b) % 10  # Simple addition modulo 10
    next_b = (a * b) % 10  # Product modulo 10
    return (next_a, next_b)


# Initial bit pair from Pi (1, 4)
initial_pair = (1, 4)

# Generate a sequence of 10 steps
sequence = recursive_growth(initial_pair, 10, example_mapping)

# Map the sequence to ASCII characters or other representations
ascii_sequence = "".join(chr(48 + a) + chr(48 + b) for a, b in sequence)

print("Generated Sequence:", sequence)
print("Mapped ASCII Sequence:", ascii_sequence)
