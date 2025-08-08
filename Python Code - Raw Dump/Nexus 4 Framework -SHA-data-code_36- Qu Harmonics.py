import hashlib
import numpy as np
import matplotlib.pyplot as plt
import statistics

# ------------------------------------------------------------------------------
# 1) Prepare SHA-256 round constants K0…K63 (first 32 bits of cube roots of primes)
# ------------------------------------------------------------------------------
PRIMES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47, 53,
    59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127, 131,
    137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223,
    227, 229, 233, 239, 241, 251, 257, 263,
    269, 271, 277, 281, 283, 293, 307, 311
]

def cube_root_fraction_bits(p):
    """Compute first 32 bits of fractional part of cube root of prime p."""
    frac = (p ** (1/3)) % 1
    return int(frac * (2**32))  # 32-bit fraction

K = [cube_root_fraction_bits(p) for p in PRIMES]  # K0…K63
Kmod256 = [k & 0xFF for k in K]

# ----------------------------------------------------------------------
# 2) Define the 13-cycle states from your feedback simulation
# ----------------------------------------------------------------------
cycle_states = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

# ------------------------------------------------------------
# A) Nearest-distance (mod 256) to each K_j for each cycle state
# ------------------------------------------------------------
print("=== A) Mod-256 Nearest-Distance to Round Constants ===")
for s in cycle_states:
    dists = [min((s - kj) % 256, (kj - s) % 256) for kj in Kmod256]
    mind = min(dists)
    js = [j for j, d in enumerate(dists) if d == mind]
    print(f"State {s:3d}: min |s−K_j| mod256 = {mind:3d}, closest j = {js}")

# ------------------------------------------------------------
# B) Repeat nearest-distance for mod 64, mod 128, mod 256
# ------------------------------------------------------------
print("\n=== B) Nearest-Distance at Various Moduli ===")
for mod in (64, 128, 256):
    print(f"\n-- Modulus = {mod} --")
    Kmod = [k % mod for k in Kmod256]
    for s in cycle_states:
        dists = [min((s - kj) % mod, (kj - s) % mod) for kj in Kmod]
        mind = min(dists)
        js = [j for j, d in enumerate(dists) if d == mind]
        print(f"State {s:3d}: min |s−K_j| mod{mod} = {mind:3d}, closest j = {js}")

# -------------------------------------------------------------------
# C) Heatmap of full 256-bit hash ints versus the 32-bit round constants
# -------------------------------------------------------------------
print("\n=== C) Heatmap of |SHA(s) - K_j| across cycle states ===")
# Compute full SHA(s) as 256-bit integers
hash_ints = []
for s in cycle_states:
    digest = hashlib.sha256(bytes([s])).digest()
    h_int = int.from_bytes(digest, 'big')
    hash_ints.append(h_int)

# Build distance matrix (states × constants)
dmat = np.zeros((len(cycle_states), len(K)), dtype=np.int64)
for i, h in enumerate(hash_ints):
    for j, kj in enumerate(K):
        # Use absolute difference, no overflow possible with Python's arbitrary precision ints
        dmat[i, j] = abs(h - kj)

# Plot
plt.figure(figsize=(10, 6))
plt.imshow(dmat, aspect='auto', cmap='viridis')
plt.colorbar(label="|SHA256(s) – K_j|")
plt.yticks(range(len(cycle_states)), cycle_states)
plt.xticks(range(0, len(K), 4), labels=range(0, len(K), 4))
plt.xlabel("Round-constant index j")
plt.ylabel("Cycle state s")
plt.title("Heatmap: Distance from SHA256(s) to Each K_j")
plt.tight_layout()
plt.show()
