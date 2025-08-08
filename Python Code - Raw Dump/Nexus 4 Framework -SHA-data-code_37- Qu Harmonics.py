import hashlib
import matplotlib.pyplot as plt

# --- Prepare constants K as before ---
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19,
          23, 29, 31, 37, 41, 43, 47, 53,
          59, 61, 67, 71, 73, 79, 83, 89,
          97, 101, 103, 107, 109, 113, 127, 131,
          137, 139, 149, 151, 157, 163, 167, 173,
          179, 181, 191, 193, 197, 199, 211, 223,
          227, 229, 233, 239, 241, 251, 257, 263,
          269, 271, 277, 281, 283, 293, 307, 311]

def cube_root_fraction_bits(p):
    frac = (p ** (1/3)) % 1
    return int(frac * (2**32))

K = [cube_root_fraction_bits(p) for p in PRIMES]

# --- Your 13-cycle states ---
cycle_states = [21, 8, 154, 65, 253, 155, 94, 203, 29, 205, 209, 153, 88]

# --- Compute full 256-bit SHA(s) ints ---
hash_ints = []
for s in cycle_states:
    h = int.from_bytes(hashlib.sha256(bytes([s])).digest(), 'big')
    hash_ints.append(h)

# --- Build a pure-Python list of lists of abs(h - K_j) ---
dist_matrix = []
for h in hash_ints:
    row = [abs(h - kj) for kj in K]
    dist_matrix.append(row)

# --- Plot with Matplotlib directly from the list of lists ---
plt.figure(figsize=(10, 6))
plt.imshow(dist_matrix, aspect='auto', cmap='viridis')
plt.colorbar(label="|SHA256(s) – K_j|")
plt.yticks(range(len(cycle_states)), cycle_states)
plt.xticks(range(0, len(K), 4), labels=range(0, len(K), 4))
plt.xlabel("Round-constant index j")
plt.ylabel("Cycle state s")
plt.title("Heatmap: Distance from SHA256(s) to Each K_j")
plt.tight_layout()
plt.show()
