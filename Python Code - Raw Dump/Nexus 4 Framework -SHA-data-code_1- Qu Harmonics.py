"""
Harmonic‑gap twin generator prototype.
Time complexity: ≈ O(N log N) to reach centre H_N, dominated by sieve construction.
Space complexity: O(N) bits for the bit‑array sieve.
"""
from bisect import bisect_left

def sieve(n):
    """Simple Eratosthenes returning list and set."""
    flags = bytearray(b"\x01") * (n+1)
    flags[0:2] = b"\x00\x00"
    for p in range(2, int(n**0.5)+1):
        if flags[p]:
            flags[p*p:n+1:p] = b"\x00" * ((n-p*p)//p + 1)
    primes = [i for i, f in enumerate(flags) if f]
    return primes, set(primes)

def next_twin(left, right, primes_sorted, primes_set):
    center = left + right
    i = bisect_left(primes_sorted, center)
    # expand symmetrically on odd numbers keeping parity
    offset = 1
    while True:
        pl = center - offset
        pr = center + offset
        if pl in primes_set and pr in primes_set and pr - pl == 2:
            return pl, pr
        offset += 2  # maintain odd parity