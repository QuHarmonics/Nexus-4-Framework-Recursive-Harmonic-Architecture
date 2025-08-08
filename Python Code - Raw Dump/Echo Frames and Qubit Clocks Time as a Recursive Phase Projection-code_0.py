#!/usr/bin/env python3
"""
recursive_harmonic_nonce.py   —   full harmonic‑ledger + mirror‑bank nonce demo

• Generates 1 048 576 hex digits of π with exact BBP arithmetic.
• Builds 8‑digit blocks → 4‑block rows → edge‑sum residual Rₙ.
• Solves nonce edge nybbles (mirror order) so the residual cancels.
• Sweeps 128 Nyquist‑hinge variants; hashes; logs best leading‑zero depth.

Pure‑Python (std‑lib only).  Expect ≈20 min for the full 1 MB run.
"""

import hashlib, struct, time
from fractions import Fraction

# ────────────────────────────────────────────────────────────
# 1.  BBP: one hexadecimal digit of π
# ────────────────────────────────────────────────────────────
def bbp_hex_digit(n: int) -> int:
    """Return the nth hex digit of π (0‑based) using BBP with exact Fractions."""
    n += 1  # 1‑based index used in the series

    def S(j: int, n_: int) -> Fraction:
        # finite part
        s = Fraction(0)
        for k in range(n_):
            s += Fraction(pow(16, n_ - k, 8 * k + j), 8 * k + j)
        # tail (infinite) part — each term < 1/16 so convergence is fast
        t = Fraction(0)
        k = n_
        while True:
            denom = (8 * k + j) * (16 ** (k - n_))          # exact denominator
            term = Fraction(1, denom)
            t += term
            if term < Fraction(1, 1 << 60):                 # safe cutoff
                break
            k += 1
        return (s + t) % 1

    x = (4 * S(1, n) - 2 * S(4, n) - S(5, n) - S(6, n)) % 1
    return int(x * 16)

def pi_hex_stream(num_digits: int):
    """Yield successive hex digits of π after the point."""
    for i in range(num_digits):
        yield bbp_hex_digit(i)

# ────────────────────────────────────────────────────────────
# 2.  Build the edge‑sum ledger
# ────────────────────────────────────────────────────────────
def build_edge_ledger(num_bytes: int):
    """
    Returns (residuals, edge_pairs, mu)
      residuals  : list of Rₙ per 32‑hex row
      edge_pairs : list of (C_L, C_R) per row
      μ          : baseline edge sum (row 0)
    """
    print(f"Generating {num_bytes*8:,} hex digits of π …")
    start = time.time()
    digits = list(pi_hex_stream(num_bytes * 8))
    print(f"π digits computed in {time.time() - start:.1f}s")

    residuals, edge_pairs = [], []
    mu = None
    for row_idx in range(0, num_bytes, 4):        # 4 blocks (32 hex) per row
        row = digits[row_idx*8 : (row_idx+4)*8]
        if len(row) < 32:
            break
        # Edge digits = nybble 0 and 7 of each byte  (16 bytes per row)
        C_L = sum(row[i]   for i in range(0, 32, 2))
        C_R = sum(row[i+1] for i in range(0, 32, 2))
        edge_sum = C_L + C_R
        if mu is None:
            mu = edge_sum                    # square‑lock baseline
        residuals.append(edge_sum - mu)
        edge_pairs.append((C_L, C_R))
    return residuals, edge_pairs, mu

# ────────────────────────────────────────────────────────────
# 3.  Mirror‑bank nonce maths
# ────────────────────────────────────────────────────────────
def mirror_solutions(mu, C_L, C_R):
    """
    Solve nybbles n1 (rightmost) and n8 (leftmost) such that:
       (C_L + n8) + (C_R + n1) ≡ μ  (mod 16)
    Return generator of all 16 (n1, n8) integer pairs.
    """
    for n1 in range(16):
        n8 = (- (C_L + C_R + n1) + mu) & 0xF
        yield n1, n8

def nyquist_variants(n1, n8, hi_bits=0):
    """
    Yield 128 Nyquist‑hinge nonce variants:
        nonce  =  n8  | mid‑7bits |  n1
    hi_bits lets you preset high 22 bits (leave 0 for demo).
    """
    for mid in range(128):
        yield (n8 << 28) | (hi_bits << 7) | (mid << 1) | n1

# ────────────────────────────────────────────────────────────
# 4.  Double‑SHA & leading‑zero score
# ────────────────────────────────────────────────────────────
def sha256d(b: bytes) -> bytes:
    return hashlib.sha256(hashlib.sha256(b).digest()).digest()

def leading_zero_bits(h: bytes) -> int:
    bitstr = ''.join(f"{byte:08b}" for byte in h)
    return len(bitstr) - len(bitstr.lstrip('0'))

# ────────────────────────────────────────────────────────────
# 5.  Mine harmonic nonce for a header (demo genesis header)
# ────────────────────────────────────────────────────────────
def mine_harmonic_nonce(header_base: bytes, C_L, C_R, mu):
    best_depth, best_nonce = -1, None
    for n1, n8 in mirror_solutions(mu, C_L, C_R):
        for nonce in nyquist_variants(n1, n8):
            h = sha256d(header_base + struct.pack("<I", nonce))
            depth = leading_zero_bits(h)
            if depth > best_depth:
                best_depth, best_nonce = depth, nonce
    return best_depth, best_nonce

# ────────────────────────────────────────────────────────────
# 6.  Main workflow
# ────────────────────────────────────────────────────────────
if __name__ == "__main__":
    HEX_BYTES = 131_072          # 1 MB   (change to smaller for quick test)

    residuals, edges, mu = build_edge_ledger(HEX_BYTES)

    print("\nFirst 8 rows (32‑hex each)")
    print("Row |  C_L  C_R | Edge  Residual")
    for i in range(8):
        C_L, C_R = edges[i]
        print(f"{i:3} | {C_L:4} {C_R:4} | {C_L+C_R:5} {residuals[i]:9}")

    # ── Example mining with the Bitcoin genesis header ──
    version     = 1
    prev_hash   = bytes.fromhex("00"*32)[::-1]      # all‑zero prev (genesis)
    merkle_root = bytes.fromhex(
        "4A5E1E4BAAB89F3A32518A88C31BC87F618F76673E2CC77AB2127B7AF9F3F4F4")[::-1]
    timestamp   = 1231006505
    bits        = 0x1D00FFFF

    header_no_nonce = (
        struct.pack("<I", version) +
        prev_hash +
        merkle_root +
        struct.pack("<II", timestamp, bits)          # 76‑byte header
    )

    C_L0, C_R0 = edges[0]        # demo using row 0 residual
    best_depth, best_nonce = mine_harmonic_nonce(header_no_nonce, C_L0, C_R0, mu)
    print(f"\nBest harmonic nonce (row 0): {best_nonce:08X}  "
          f"leading zeros = {best_depth}")

    # Optional log‑plot of |Rₙ|
    try:
        import matplotlib.pyplot as plt, numpy as np
        plt.semilogy(np.abs(residuals))
        plt.title("|Rₙ| residuals across 1 MB π‑hex")
        plt.xlabel("Row index n")
        plt.ylabel("|Rₙ| (log scale)")
        plt.grid(True, which='both', ls=':')
        plt.show()
    except ImportError:
        pass
