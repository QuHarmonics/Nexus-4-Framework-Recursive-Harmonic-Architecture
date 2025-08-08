import math

def is_prime(n):
    if n < 2: return False
    if n in (2, 3): return True
    if n % 2 == 0 or n % 3 == 0: return False
    for i in range(5, int(math.isqrt(n)) + 1, 6):
        if n % i == 0 or n % (i + 2) == 0:
            return False
    return True

def bbp_delta(n, k_max=4):
    return int(sum((16 ** (1 - k)) / (8 * k + n % 7 + 1) for k in range(1, k_max + 1))) + 1

def generate_bbp_twin_primes(limit, filepath):
    twin_primes = []
    with open(filepath, 'w') as f:
        current = 3
        while current < limit:
            if is_prime(current) and is_prime(current + 2):
                pair = (current, current + 2)
                twin_primes.append(pair)
                f.write(f"{pair}\n")
            current += bbp_delta(current)
    return filepath, len(twin_primes)

    # Set the limit and file path
limit = 10_000_000
output_path = "D://twin_primes_output.txt"

# Run and write to file
generate_bbp_twin_primes(limit, output_path)