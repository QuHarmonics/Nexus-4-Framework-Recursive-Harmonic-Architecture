import hashlib
import math
import numpy as np
import matplotlib.pyplot as plt

# --- Prepare the SHA-256 round constants K_j (first 32 bits of cube roots of primes) ---
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
    frac = (p ** (1/3)) % 1
    return int(frac * (2**32))

K = [cube_root_fraction_bits(p) for p in PRIMES]

# --- Your 13-cycle states (from the Δ-phase loop) ---
cycle_states = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

# --- Compute the full 256-bit SHA256 digest ints for each state ---
hash_ints = []
for s in cycle_states:
    digest = hashlib.sha256(bytes([s])).digest()
    hash_ints.append(int.from_bytes(digest, 'big'))

# --- Build a float64 matrix of log10(dist+1) to avoid overflow ---
n_states = len(cycle_states)
n_consts = len(K)
log_dmat = np.zeros((n_states, n_consts), dtype=np.float64)

for i, h in enumerate(hash_ints):
    for j, kj in enumerate(K):
        d = abs(h - kj)
        log_dmat[i, j] = math.log10(d + 1)  # +1 so that identical values map to log10(1)=0

# --- Plot the heatmap ---
plt.figure(figsize=(10, 6))
im = plt.imshow(log_dmat, aspect='auto', cmap='magma')
plt.colorbar(im, label="log10(|SHA256(s) – K_j| + 1)")
plt.yticks(range(n_states), cycle_states)
plt.xticks(range(0, n_consts, 4), labels=range(0, n_consts, 4))
plt.xlabel("Round‐constant index j")
plt.ylabel("Cycle state s")
plt.title("Log-scaled Heatmap of |SHA256(s) – K_j| Across 13-Cycle States")
plt.tight_layout()
plt.show()
