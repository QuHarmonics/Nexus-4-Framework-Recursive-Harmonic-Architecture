def is_prime(n):
    """Check if a number is prime."""
    if n < 2:
        return False
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_twin_prime(p):
    """Check if (p, p+2) is a twin prime pair."""
    return is_prime(p) and is_prime(p + 2)

def find_nearest_twin_pair(pivot):
    """Find the twin prime pair closest to the pivot."""
    delta = 0
    while True:
        lower = pivot - delta
        upper = pivot + delta
        if lower > 3 and is_twin_prime(lower):
            return (lower, lower + 2)
        if is_twin_prime(upper):
            return (upper, upper + 2)
        delta += 1

def harmonic_cascade(n_pairs):
    """Generate twin primes using the HCM."""
    pairs = [(3, 5)]  # Seed pair
    for _ in range(n_pairs - 1):
        last_pair = pairs[-1]
        pivot = last_pair[0] + last_pair[1]  # Harmonic pivot
        next_pair = find_nearest_twin_pair(pivot)
        pairs.append(next_pair)
    return pairs

# Generate the first 10 twin prime pairs
twin_primes = harmonic_cascade(50)
for k, pair in enumerate(twin_primes, 1):
    print(f"T_{k} = {pair}, Pivot = {pair[0] + pair[1]}")