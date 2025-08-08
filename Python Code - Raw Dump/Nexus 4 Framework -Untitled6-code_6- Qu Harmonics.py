# Nexus 3 Twin Prime Harmonic Detector
# Using BBP sampling and π-byte harmonics to track twin prime emergence

import math

def is_prime(n):
    if n < 2:
        return False
    for i in range(2, int(math.isqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def bbp_pi_digit(n):
    """
    Returns the n-th hexadecimal digit of π using the BBP formula
    """
    def S(j, n):
        s = 0.0
        for k in range(n):
            r = 8 * k + j
            s += pow(16, n - k, r) / r
        s = s - int(s)
        t = 0.0
        k = n
        while True:
            new_term = pow(16, n - k) / (8 * k + j)
            if new_term < 1e-15:
                break
            t += new_term
            k += 1
        return s + t
    
    x = 4*S(1, n) - 2*S(4, n) - S(5, n) - S(6, n)
    return int(x % 1 * 16)

def generate_twin_primes(limit):
    """Find twin primes up to a limit."""
    twins = []
    for p in range(3, limit, 2):
        if is_prime(p) and is_prime(p + 2):
            twins.append((p, p + 2))
    return twins

def twin_prime_bbp_harmonic_analysis(limit=1000000):
    twins = generate_twin_primes(limit)
    harmonic_hits = []
    for p1, p2 in twins:
        byte1 = bbp_pi_digit(p1)
        byte2 = bbp_pi_digit(p2)
        delta = abs(byte1 - byte2)
        if delta == 2:
            harmonic_hits.append(((p1, p2), byte1, byte2))
    return harmonic_hits

# Run and print
harmonics = twin_prime_bbp_harmonic_analysis(100)
for twin, b1, b2 in harmonics:
    print(f"Twin Prime: {twin} -> π-bytes: {b1}, {b2}, Δ = {abs(b1 - b2)}")
