import numpy as np
import matplotlib.pyplot as plt
import hashlib

# === SETTINGS ===
input_data = b"hello world"           # Input message
modulus = 64                          # 6-bit harmonic register (mod 64)
samples = 256                         # Number of sequential states

# === HASH-BASED RESIDUE FUNCTION ===
def sha_residue(state, input_bytes, mod=64):
    """Compute SHA-256 of (state || input), return residue mod M"""
    msg = bytes([state % 256]) + input_bytes
    h = hashlib.sha256(msg).digest()
    h_int = int.from_bytes(h, 'big')
    return h_int % mod

# === BUILD RESIDUE SEQUENCE ===
residues = np.array([sha_residue(i, input_data, mod=modulus) for i in range(samples)])

# === BUILD HEXICON GRID (MOD 64 MAPPED TO 8×8 HEXAGONAL TOPOLOGY) ===
grid_size = 8
grid = np.zeros((grid_size, grid_size))

for i, r in enumerate(residues):
    x = r % grid_size
    y = r // grid_size
    grid[y, x] += 1  # Increment density

# === PLOT HEXICON RESIDUE GRID ===
plt.figure(figsize=(8, 7))
plt.imshow(grid, cmap='viridis', interpolation='nearest')
plt.colorbar(label='Residue Frequency')
plt.title("Hexicon Grid — SHA Residues in 6-Bit Harmonic Space", fontsize=14)
plt.xlabel("Residue Column (mod 8)")
plt.ylabel("Residue Row (mod 8)")
plt.xticks(np.arange(grid_size))
plt.yticks(np.arange(grid_size))
plt.grid(False)
plt.tight_layout()
plt.show()
