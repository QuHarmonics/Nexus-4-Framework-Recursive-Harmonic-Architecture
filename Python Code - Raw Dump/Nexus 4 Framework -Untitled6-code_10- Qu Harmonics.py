import math

# Prime check
def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

# Binary length
def binary_length(n):
    return len(bin(n)[2:])

# BBP-like inspired step function
def bbp_delta(n):
    step = 0
    for k in range(1, 5):  # Small bounded sum for simplicity
        step += (16 ** (1 - k)) / (8 * k + n % 7 + 1)
    return int(step) + 1  # Ensure forward movement

# BBP-modulated twin prime generator with file output
def generate_bbp_twin_primes_to_file(limit, filepath):
    twin_primes = []
    current = 3
    with open(filepath, 'w') as f:
        while current < limit:
            if is_prime(current) and is_prime(current + 2):
                twin_pair = (current, current + 2)
                twin_primes.append(twin_pair)
                f.write(f"{twin_pair}\n")
            step = bbp_delta(current)
            current += step
    return filepath, len(twin_primes)

# Set the limit and file path
limit = 10_000_000
output_path = "D://twin_primes_output.txt"

# Run and write to file
generate_bbp_twin_primes_to_file(limit, output_path)
