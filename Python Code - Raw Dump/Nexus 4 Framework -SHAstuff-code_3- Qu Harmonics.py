import hashlib
import re
from collections import defaultdict

def sha256_hexdigest(s: str) -> str:
    return hashlib.sha256(s.encode('utf-8')).hexdigest()

def ascii_hexify(hex_str: str) -> str:
    """Convert each hex digit into its two‑char ASCII hex code."""
    return ''.join(f"{ord(c):02x}" for c in hex_str)

def reverse_nibbles(ascii_hex: str) -> str:
    """Reverse the string one hex digit (nibble) at a time."""
    return ascii_hex[::-1]

def ascii_dehexify(rev: str) -> str:
    """Turn reversed ASCII‑hex back into characters (for debugging)."""
    return ''.join(chr(int(rev[i:i+2], 16)) 
                   for i in range(0, len(rev), 2))

def trailing_zero_count(n: int) -> int:
    """Count how many times 16 divides the integer (i.e. trailing hex zeros)."""
    count = 0
    while n and n % 16 == 0:
        n //= 16
        count += 1
    return count

def is_power_of_16(n: int) -> bool:
    return n > 0 and (n & (n - 1) == 0) and all(((n.bit_length() - 1) % 4) == 0 for _ in [0])

def decaying_wave_score(deltas: list[int]) -> float:
    """
    A simple “decay” metric: how quickly the sequence of absolute differences
    drops off. We'll just take the ratio between the first and last few.
    """
    if len(deltas) < 4:
        return 0.0
    head = sum(abs(x) for x in deltas[:2]) / 2
    tail = sum(abs(x) for x in deltas[-2:]) / 2
    return (head - tail) / (head + 1e-12)

class ResonanceDifferentiator:
    def __init__(self, target: str):
        # Accept target as hex or raw string
        self.target_hash = (target if re.fullmatch(r"[0-9a-f]{64}", target)
                            else sha256_hexdigest(target))
        # Precompute reversed ASCII‑hex for target
        t_ascii = ascii_hexify(self.target_hash)
        self.target_rev = reverse_nibbles(t_ascii)
        self.target_int = int(self.target_rev, 16)

    def analyze(self, guess: str) -> dict:
        g_hash = sha256_hexdigest(guess)
        g_ascii = ascii_hexify(g_hash)
        g_rev = reverse_nibbles(g_ascii)
        g_int = int(g_rev, 16)

        delta = abs(self.target_int - g_int)
        # Compute a rough “wave” of differences (for decay scoring)
        # Here we break the reversed strings into word‑sized chunks:
        chunk_size = 8
        deltas = []
        for i in range(0, len(self.target_rev), chunk_size):
            t_chunk = int(self.target_rev[i:i+chunk_size], 16)
            g_chunk = int(g_rev[i:i+chunk_size], 16)
            deltas.append(t_chunk - g_chunk)

        return {
            "guess": guess,
            "hash": g_hash,
            "delta": delta,
            "trailing_zeros": trailing_zero_count(delta),
            "power_of_16": is_power_of_16(delta),
            "decay_score": decaying_wave_score(deltas),
        }

    def rank(self, guesses: list[str]) -> list[dict]:
        results = [self.analyze(g) for g in guesses]
        # Sort by best harmonicity: more zeros, power‑of‑16 first, then decay_score
        return sorted(
            results,
            key=lambda r: (
                -r["power_of_16"],       # True first
                -r["trailing_zeros"],    # More zeros
                -r["decay_score"],       # Stronger decay
            )
        )

if __name__ == "__main__":
    target = "The ultimate truth."
    engine = ResonanceDifferentiator(target)

    candidates = [
        "hello", "Hello", "hallo", "hullo",
        "0xdeadbeef", "42", "resonance", "resonator"
    ]
    ranked = engine.rank(candidates)

    from pprint import pprint
    pprint(ranked, width=120)
