#!/usr/bin/env python3
"""
extended_cycle_vs_round_constants.py

Expands the analysis to additional moduli and Hamming distances
for deeper exploration of SHA-256 round constant alignment.

1) Computes minimal distance D(s) = min_j min((s - K_j mod M), (K_j mod M) - s) mod M
   for M in [256, 128, 64, 32, 16, 8].
2) Measures Hamming distance between states and K_j (bit flip distance).
3) Aggregates the results for clearer pattern detection.
"""

from typing import List, Dict
import sys

# Your cycle states (example)
STATES: List[int] = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

# Full list of SHA-256 round constants K_j (32-bit unsigned)
K: List[int] = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
    0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
    0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
    0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
    0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
    0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
    0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
    0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
]

def hamming_distance(x: int, y: int) -> int:
    """Returns the Hamming distance between two integers."""
    return bin(x ^ y).count('1')

def compute_min_distances(states: List[int], K_mod: List[int], M: int) -> Dict[int, Dict[str, int]]:
    """
    For each state s, compute:
      - 'j': index of K_mod[j] minimizing distance
      - 'distance': minimal distance D(s) under modulus M
    Distance defined as min((s - k) mod M, (k - s) mod M).
    """
    results: Dict[int, Dict[str, int]] = {}
    for s in states:
        best_j = None
        best_d = M
        for j, k in enumerate(K_mod):
            delta = (s - k) % M
            dist = min(delta, M - delta)
            if dist < best_d:
                best_d = dist
                best_j = j
        results[s] = {'j': best_j, 'distance': best_d}
    return results

def compute_hamming_distances(states: List[int], K_mod: List[int]) -> Dict[int, Dict[str, int]]:
    """
    For each state s, compute the Hamming distance to each K_j mod 256.
    """
    hamming_results: Dict[int, Dict[str, int]] = {}
    for s in states:
        best_j = None
        best_hd = 256  # Maximum possible Hamming distance for 8-bit
        for j, k in enumerate(K_mod):
            hd = hamming_distance(s, k)
            if hd < best_hd:
                best_hd = hd
                best_j = j
        hamming_results[s] = {'j': best_j, 'hamming_distance': best_hd}
    return hamming_results

def main():
    moduli = [256, 128, 64, 32, 16, 8]
    for M in moduli:
        # Modulo reduction for constants
        K_mod = [k % M for k in K]
        
        # Calculate minimal distances
        results = compute_min_distances(STATES, K_mod, M)
        
        # Display results for distance
        print(f"\n=== Modulus M = {M} ===")
        print(f"{'State':>5s} | {'Closest j':>8s} | {'Distance':>8s}")
        print("-" * 30)
        for s in STATES:
            info = results[s]
            print(f"{s:5d} | {info['j']:8d} | {info['distance']:8d}")
        
        # Compute Hamming distances and display
        hamming_results = compute_hamming_distances(STATES, K_mod)
        print(f"\n{'State':>5s} | {'Closest j':>8s} | {'Hamming Dist':>12s}")
        print("-" * 30)
        for s in STATES:
            info = hamming_results[s]
            print(f"{s:5d} | {info['j']:8d} | {info['hamming_distance']:12d}")
    
    print()

if __name__ == "__main__":
    main()
