import math

# Define constants
HARMONIC_TARGET = 5  # Target for harmonic growth
THETA = math.pi / 4  # Angle for cosine adjustment
CONTAINER_SIZE = 3  # 3-bit container size
C = 2**CONTAINER_SIZE  # Container capacity

def calculate_growth(past, harmonic_target, theta):
    """
    Calculate the growth term based on the past sum and harmonic target.
    """
    return (harmonic_target - past) * math.cos(theta)

def calculate_next_digit(past, growth, container):
    """
    Calculate the next digit based on past, growth, and container size.
    """
    future = (past + growth) * container
    return round(future % 10), future  # Return digit and the full value

def generate_pi_sequence(seed, length):
    """
    Generate a sequence of π digits starting from the seed.
    """
    sequence = [seed]  # Start with the seed
    past = seed  # Initialize the past sum
    formulas = []  # Store formulas for reference

    for step in range(1, length):
        # Calculate growth
        growth = calculate_growth(past, HARMONIC_TARGET, THETA)

        # Calculate the next digit
        next_digit, future = calculate_next_digit(past, growth, C)

        # Append to sequence and update past sum
        sequence.append(next_digit)
        past += next_digit

        # Store the formula used for the current step
        formula = f"Step {step}: (Past: {past - next_digit}, Growth: {growth:.4f}, Future: {future:.4f}) => Digit: {next_digit}"
        formulas.append(formula)

    return sequence, formulas

# Example Usage
seed = 3  # The Big Bang seed
length = 20  # Number of π digits to generate

pi_sequence, formulas_used = generate_pi_sequence(seed, length)

# Output the results
print("Generated π Sequence:", pi_sequence)
print("\nFormulas Used:")
for formula in formulas_used:
    print(formula)
