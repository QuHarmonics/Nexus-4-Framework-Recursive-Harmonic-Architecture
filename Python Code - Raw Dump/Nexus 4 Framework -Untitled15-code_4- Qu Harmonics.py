import math

def calculate_delta(a, b):
    """
    Calculate the delta (difference between current and past values).
    """
    return b - a

def calculate_container_size(a, b, delta):
    """
    Determine the container size based on the binary length of (A + B + Delta).
    """
    combined_state = a + b + delta
    binary_length = len(bin(combined_state)) - 2  # Exclude the '0b' prefix
    return 2 ** binary_length

def calculate_future_value(a, b, delta, container_size):
    """
    Calculate the future value (F) based on A, B, Delta, and container size.
    """
    return (a + b + delta) * container_size

def calculate_next_digit(future_value):
    """
    Determine the next digit based on the binary length of the future value.
    """
    binary_representation = bin(int(abs(future_value)))
    return len(binary_representation) - 2  # Exclude the '0b' prefix

def generate_sequence(a, b, iterations):
    """
    Generate a sequence starting with A and B using the recursive Nexus formula.
    """
    sequence = [a, b]  # Start with A and B

    for _ in range(iterations - 2):  # Two initial values already in the sequence
        # Step 1: Calculate Delta
        delta = calculate_delta(sequence[-2], sequence[-1])

        # Step 2: Determine Container Size
        container_size = calculate_container_size(sequence[-2], sequence[-1], delta)

        # Step 3: Compute Future Value
        future_value = calculate_future_value(sequence[-2], sequence[-1], delta, container_size)

        # Step 4: Derive Next Digit
        next_digit = calculate_next_digit(future_value)

        # Update sequence
        sequence.append(next_digit)

    return sequence

# Parameters
a = 1  # Past value
b = 4  # Current value
iterations = 10  # Number of digits to generate

# Generate the sequence
sequence = generate_sequence(a, b, iterations)

# Output the results
print("Generated Sequence:", sequence)
