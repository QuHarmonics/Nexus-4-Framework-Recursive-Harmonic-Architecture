#!/usr/bin/env python3
"""
twin_prime_analysis.py

Identify primes and twin primes within a list of nibble sums.

1) Checks each value for primality.
2) Finds all twin-prime pairs (pairs of primes differing by 2).
3) Reports primes and twin-prime pairs.

Usage:
    python twin_prime_analysis.py
"""

from typing import List, Tuple

# The nibble sums you provided
NIBBLE_SUMS: List[int] = [5, 6, 11, 15, 16, 12, 15, 17, 19, 21, 11, 13, 15]

def is_prime(n: int) -> bool:
    """Return True if n is prime, False otherwise."""
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def find_primes(values: List[int]) -> List[int]:
    """Return a list of the primes in values (in original order, duplicates preserved)."""
    return [v for v in values if is_prime(v)]

def find_twin_primes(primes: List[int]) -> List[Tuple[int,int]]:
    """
    Given a list of primes, return all twin-prime pairs (p, p+2)
    that both appear in the list.
    Each pair is reported once, with p < p+2.
    """
    prime_set = set(primes)
    twins = []
    for p in primes:
        if p + 2 in prime_set:
            twins.append((p, p + 2))
    # Remove duplicates (e.g., both (11,13) and (13,11) scenarios)
    unique_twins = []
    seen = set()
    for a, b in twins:
        if (a, b) not in seen and (b, a) not in seen:
            unique_twins.append((a, b))
            seen.add((a, b))
    return unique_twins

def main():
    print("Nibble sums:", NIBBLE_SUMS)
    primes = find_primes(NIBBLE_SUMS)
    print("Primes in nibble sums:", primes)

    twins = find_twin_primes(primes)
    if twins:
        print("Twin-prime pairs found:")
        for a, b in twins:
            print(f"  ({a}, {b})")
    else:
        print("No twin-prime pairs found.")

if __name__ == "__main__":
    main()
