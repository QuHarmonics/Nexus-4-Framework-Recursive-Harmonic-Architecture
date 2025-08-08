import math

# --- BBP Delta Step Function ---
def bbp_delta(n, kmax=4):
    step = sum(16**(1 - k) / (8 * k + n % 7 + 1) for k in range(1, kmax + 1))
    return int(math.floor(step)) + 1

# --- Twin Prime Generator Using BBP Step ---
def twin_primes_bbp(limit, kmax=4):
    pairs = []
    n = 3
    while n < limit:
        if is_prime(n) and is_prime(n + 2):
            pairs.append((n, n + 2))
        n += bbp_delta(n, kmax)
    return pairs

# --- Primality Check ---
def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(math.isqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True

# --- Export Twin Primes to File ---
def export_twin_primes(limit, filepath):
    data = twin_primes_bbp(limit)
    with open(filepath, 'w') as f:
        for a, b in data:
            f.write(f"{a}\t{b}\n")
    return {"FilePath": filepath, "Count": len(data)}

# --- Usage Example ---
# Uncomment the following line to run:
# export_twin_primes(10**9, "d:/twin_primes_output.txt")
