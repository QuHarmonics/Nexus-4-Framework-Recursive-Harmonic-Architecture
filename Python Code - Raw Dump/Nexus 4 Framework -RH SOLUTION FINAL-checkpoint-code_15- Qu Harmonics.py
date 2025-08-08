import mpmath

# Step 1: Calculate Pi to 1000 decimal places
mpmath.mp.dps = 1000  # Set precision
pi_digits = str(mpmath.pi)[2:]  # Get π digits after "3."

# Step 2: Function to generate sequence with padding
def generate_sequence_with_padding(pi_digits):
    sequence = []
    for i in range(1, len(pi_digits)):
        a = int(pi_digits[i - 1])  # Current number
        b = int(pi_digits[i])      # Previous number
        c_padding = len(str(a**2 + b**2))  # Length of a^2 + b^2
        
        # Add current number + padding zeros
        padded_entry = f"{a}" + "0" * c_padding
        sequence.append(padded_entry)
    sequence.append(pi_digits[-1])  # Add the last digit of π
    return "".join(sequence)

# Step 3: Generate and print the padded sequence
padded_sequence = generate_sequence_with_padding(pi_digits)
print(padded_sequence)

