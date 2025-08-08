#!/usr/bin/env python3
"""
enhanced_twin_prime_analysis.py

Identify primes, twin primes, and solitary primes within a list of nibble sums.

1) Checks each value for primality.
2) Finds all twin-prime pairs (pairs of primes differing by 2).
3) Reports:
   • All primes in the list (with duplicates).
   • Twin-prime pairs found.
   • Solitary primes (primes that are not in any twin pair).

Usage:
    python enhanced_twin_prime_analysis.py
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

def find_twin_primes(primes: List[int]) -> List[Tuple[int, int]]:
    """
    Given a list of primes, return all twin-prime pairs (p, p+2)
    that both appear in the list. Each pair is reported once with p < p+2.
    """
    prime_set = set(primes)
    twins = []
    for p in prime_set:
        if p + 2 in prime_set:
            twins.append((p, p + 2))
    return sorted(twins)

def find_solitary_primes(primes: List[int], twins: List[Tuple[int, int]]) -> List[int]:
    """
    Return primes that are not part of any twin-prime pair.
    Duplicates preserved in order of first appearance.
    """
    twin_members = {p for pair in twins for p in pair}
    solitary = []
    for p in primes:
        if p not in twin_members and p not in solitary:
            solitary.append(p)
    return solitary

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

    solitary = find_solitary_primes(primes, twins)
    if solitary:
        print("Solitary primes (not part of any twin):", solitary)
    else:
        print("No solitary primes.")

if __name__ == "__main__":
    main()
