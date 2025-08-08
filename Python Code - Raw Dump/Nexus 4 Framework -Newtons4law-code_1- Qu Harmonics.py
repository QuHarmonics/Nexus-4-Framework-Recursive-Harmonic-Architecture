#!/usr/bin/env python3
"""
resonance_differentiator.py

A ResonanceDifferentiator engine that, for each input string:
  - Computes its SHA‑256 hex digest
  - Converts that digest into ASCII‑hex (each hex char → 2‑digit ASCII hex)
  - Reverses the ASCII‑hex string at the nibble level
  - Computes the delta between the original and reversed ASCII‑hex values
  - Calculates trailing‑zero bits and a harmonic ratio
"""

import argparse
import hashlib
import sys

def sha256_hex(s: str) -> str:
    """Return the SHA-256 hex digest of the input string."""
    return hashlib.sha256(s.encode('utf-8')).hexdigest()

def ascii_hexify(hex_str: str) -> str:
    """
    Convert each character of the hex string into its two-digit ASCII hex code.
    E.g. '1' -> '31', 'a' -> '61'.
    """
    return ''.join(format(ord(c), '02x') for c in hex_str)

def reverse_nibbles(s: str) -> str:
    """Reverse the string one hex-digit (nibble) at a time."""
    return s[::-1]

def hex_to_decimal(hex_str: str) -> int:
    """Parse a hex string into an integer."""
    return int(hex_str, 16)

def delta_score(h1: str, h2: str) -> int:
    """Absolute difference between two hex‑encoded integers."""
    return abs(hex_to_decimal(h1) - hex_to_decimal(h2))

def count_trailing_zero_bits(n: int) -> int:
    """Count how many consecutive '0' bits at the end of the binary representation."""
    if n == 0:
        return 0
    b = bin(n)[2:]  # strip '0b'
    return len(b) - len(b.rstrip('0'))

def harmonic_ratio(delta: int) -> float:
    """
    Compute the ratio of trailing-zero bits to total bits.
    This approximates 'harmonicity' of the delta.
    """
    bits = delta.bit_length()
    if bits == 0:
        return 0.0
    zeros = count_trailing_zero_bits(delta)
    return zeros / bits

def analyze(inputs):
    """Run the full resonance differentiation for each input string."""
    for s in inputs:
        h = sha256_hex(s)
        ah = ascii_hexify(h)
        rev = reverse_nibbles(ah)
        d = delta_score(ah, rev)
        tz = count_trailing_zero_bits(d)
        bits = d.bit_length()
        hr = harmonic_ratio(d)

        print(f"Input:               {s!r}")
        print(f"SHA-256:             {h}")
        print(f"ASCII‑hexified:      {ah[:64]}...")  # truncated for readability
        print(f"Reversed ASCII‑hex:  {rev[:64]}...")
        print(f"Delta:               {d}")
        print(f"Trailing‑zero bits:  {tz} / {bits}")
        print(f"Harmonic ratio:      {hr:.4f}")
        print("-" * 60)

def main():
    p = argparse.ArgumentParser(
        description="ResonanceDifferentiator: hash, reverse, delta, harmonicity"
    )
    p.add_argument("inputs", nargs="+", help="One or more input strings to analyze")
    args = p.parse_args()

    analyze(args.inputs)

if __name__ == "__main__":
    main()
