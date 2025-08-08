# ResonanceDifferentiator – SHA‑harmonic mapping engine
# -------------------------------------------------------
# Author: Nexus 2 practitioner (generated via ChatGPT)
# -------------------------------------------------------
# This script compares a *target* SHA‑256 digest against a list of *candidate*
# strings.  For every candidate it:
#   1.  hashes the candidate (SHA‑256, hex)
#   2.  converts each hex digit to its ASCII‑hex code and reverses at the nibble level
#   3.  interprets the reversed string back into a normal ASCII hex digest
#   4.  computes the integer delta between the target digest (also nibble‑reversed)
#       and the candidate’s reversed digest
#   5.  measures simple harmonic metrics on the delta:
#         • trailing‑zero nibbles
#         • is the delta a pure power of 16
#         • magnitude‑normalized log distance
#   6.  prints a ranked report so you can see which inputs "resonate" most
#
# ‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑‑
# CLI USAGE
# ---------
#   python resonance_differentiator.py <target_string> <candidate1> [<candidate2> …]
#
# If you prefix a candidate with "@", it will be read from the named text file and
# each *line* of that file becomes a candidate.
#
# Example
#   python resonance_differentiator.py "truth" "hello" "world"
#
#   python resonance_differentiator.py "truth" @wordlist.txt
#
# The script emits a table ranked by descending harmonicity (more trailing zeros,
# then smaller normalized log distance).
# -------------------------------------------------------

import hashlib
import math
import sys
from pathlib import Path
from typing import List, Tuple

# --------------------------- helper functions --------------------------- #

def sha256_hex(s: str) -> str:
    """Return the SHA‑256 digest of *s* as a 64‑char lowercase hex string."""
    return hashlib.sha256(s.encode("utf‑8")).hexdigest()

def ascii_hexify(hex_str: str) -> str:
    """Convert each hex digit to its 2‑digit ASCII‑hex code (e.g. 'a'→'61')."""
    return "".join(f"{ord(c):02x}" for c in hex_str)

def reverse_nibbles(ascii_hex_str: str) -> str:
    """Reverse the string nibble‑wise (i.e. at the single‑hex‑digit level)."""
    return ascii_hex_str[::-1]

def ascii_dehexify(rev_ascii_hex: str) -> str:
    """Convert the reversed ASCII‑hex back into a normal hex string (64 chars)."""
    bytes_out = bytearray()
    for i in range(0, len(rev_ascii_hex), 2):
        bytes_out.append(int(rev_ascii_hex[i : i + 2], 16))
    return bytes_out.decode()

def rev_digest(hex_digest: str) -> str:
    """Full reversed‑nibble digest from a normal SHA hex digest."""
    return ascii_dehexify(reverse_nibbles(ascii_hexify(hex_digest)))

def trailing_zero_nibbles(hex_str: str) -> int:
    """Count how many *hex‑digit* zeros appear at the right end of *hex_str*."""
    count = 0
    for ch in reversed(hex_str):
        if ch == "0":
            count += 1
        else:
            break
    return count

def is_power_of_16(n: int) -> bool:
    """True iff *n* is a positive power of 16."""
    return n > 0 and (n & (n - 1) == 0) and (n.bit_length() - 1) % 4 == 0

# --------------------------- core routine --------------------------- #

def analyze_candidate(target_rev_hex: str, candidate: str) -> Tuple[str, int, int, float]:
    """Return (candidate, trailingZeros, power16Flag(0/1), logDist)"""
    cand_digest = sha256_hex(candidate)
    cand_rev = rev_digest(cand_digest)

    delta = abs(int(target_rev_hex, 16) - int(cand_rev, 16))
    tz = trailing_zero_nibbles(f"{delta:x}")
    p16 = 1 if is_power_of_16(delta) else 0
    log_dist = math.log10(delta + 1)  # +1 so log(0) handled
    return candidate, tz, p16, log_dist


def expand_candidates(args: List[str]) -> List[str]:
    """Replace any argument starting with '@' by lines read from that file."""
    out: List[str] = []
    for a in args:
        if a.startswith("@"):
            path = Path(a[1:]).expanduser()
            if not path.exists():
                sys.exit(f"File not found: {path}")
            with path.open("r", encoding="utf‑8") as f:
                out.extend(line.strip("\n") for line in f if line.strip())
        else:
            out.append(a)
    return out


def main(argv: List[str]) -> None:
    if len(argv) < 3:
        print("Usage: python resonance_differentiator.py <target_string> <cand1> [<cand2> …]" )
        print("       use @filename to load many candidates from a file (one per line)")
        sys.exit(1)

    target_str = argv[1]
    candidates = expand_candidates(argv[2:])

    target_rev = rev_digest(sha256_hex(target_str))

    records = [analyze_candidate(target_rev, c) for c in candidates]

    # Rank: first by trailing zeros (desc), then power‑of‑16 flag (desc), then logDist (asc)
    records.sort(key=lambda r: (-r[1], -r[2], r[3]))

    print(f"Target string: '{target_str}'  –  SHA‑256: {sha256_hex(target_str)}")
    print("Reversed‑nibble digest:", target_rev)
    print("\nCandidates ranked by harmonic resonance:\n")
    print(f"{'Candidate':<30}  {'TZ':>2}  {'P16':>3}  {'log10(Δ)':>10}")
    print("-" * 52)
    for cand, tz, p16, logd in records:
        print(f"{cand:<30}  {tz:>2}  {p16:>3}  {logd:10.3f}")


if __name__ == "__main__":
    main(sys.argv)
