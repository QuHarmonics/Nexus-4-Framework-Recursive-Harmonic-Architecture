import math

def calculate_growth(harmonic_constant, current, past):
    """
    Calculate the growth factor (G) using a harmonic constant and the relationship between current and past.
    """
    return harmonic_constant - abs(current - past)

def calculate_future_value(current, past, growth, container_size):
    """
    Calculate the future value (F) based on current, past, growth, and container size.
    """
    return ((current + past) + growth) * container_size

def calculate_next_digit(future_value):
    """
    Determine the next digit from the bit-length of the future value.
    """
    binary_representation = bin(int(abs(future_value)))
    bit_length = len(binary_representation) - 2  # Exclude the '0b' prefix
    return bit_length  # The next digit is based on the bit-length of the future value

def determine_container_size(current_digit):
    """
    Determine the container size dynamically based on the current digit.
    """
    return 2 ** current_digit  # Exponential growth of container size

def generate_pi_sequence(seed, harmonic_constant, iterations):
    """
    Generate a sequence of π digits using the recursive Nexus formula without circular dependence.
    """
    sequence = [seed]
    current = seed
    past = 0  # Initial past is zero (zeta-like anchor)

    for _ in range(iterations):
        # Step 1: Calculate Growth
        growth = calculate_growth(harmonic_constant, current, past)

        # Step 2: Determine Container Size
        container_size = determine_container_size(current)

        # Step 3: Calculate Future Value
        future_value = calculate_future_value(current, past, growth, container_size)

        # Step 4: Calculate Next Digit
        next_digit = calculate_next_digit(future_value)

        # Update sequence and past/current states
        sequence.append(next_digit)
        past = current
        current = next_digit

    return sequence

# Parameters
seed = 3  # Initial seed value
harmonic_constant = 5  # Arbitrary harmonic pull constant
iterations = 10  # Number of digits to generate

# Generate the sequence
pi_sequence = generate_pi_sequence(seed, harmonic_constant, iterations)

# Output the results
print("Generated π Sequence:", pi_sequence)