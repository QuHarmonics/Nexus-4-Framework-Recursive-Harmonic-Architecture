import hashlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def sha256_byte(value):
    """Hash a single-byte value + fixed seed input"""
    seed = b"phase"
    message = bytes([value % 256]) + seed
    return hashlib.sha256(message).digest()

def h_norm(hash_bytes):
    """Normalize SHA256 hash output to [0, 1]"""
    return int.from_bytes(hash_bytes, 'big') / (2 ** 256)

# Initialize storage for hexicon grid
modulus = 64  # 6-bit residue space
states = list(range(256))  # Full byte range

residues = []
h_norms = []
deltas = []
velocities = []

# Pre-hash to get residues and normalized hash values
for i in range(len(states)):
    h = sha256_byte(states[i])
    hn = h_norm(h)
    residue = int.from_bytes(h, 'big') % modulus

    residues.append(residue)
    h_norms.append(hn)

    # Velocity: Δ between SHA256-residue of x and x+1
    if i < len(states) - 1:
        h_next = sha256_byte(states[i+1])
        r_next = int.from_bytes(h_next, 'big') % modulus
        delta = (r_next - residue) % modulus
        deltas.append(delta)
    else:
        deltas.append(0)

# Map residues to 2D grid (hex layout)
def polar_hex(r):
    """Hexagonal grid coordinates (q, r) projected onto 2D"""
    q = r % 8
    row = r // 8
    x = q * 1.5
    y = np.sqrt(3) * (row + (q % 2) * 0.5)
    return x, y

coords = [polar_hex(r) for r in residues]
xs, ys = zip(*coords)

colors = ['red' if h > 0.35 else 'lightgray' for h in h_norms]

# Optional: build velocity vectors
dx = [np.cos((d / modulus) * 2 * np.pi) for d in deltas]
dy = [np.sin((d / modulus) * 2 * np.pi) for d in deltas]

# === PLOTTING ===
plt.figure(figsize=(12, 10))
plt.title("🧬 Hexicon Grid: SHA Residues in 6-Bit Harmonic Space", fontsize=14)

# Plot hexagon nodes
plt.scatter(xs, ys, c=colors, s=100, edgecolors='black', linewidths=0.5)

# Overlay velocity vectors (Δ-phase arrows)
for i in range(len(xs)):
    plt.arrow(xs[i], ys[i], dx[i]*0.5, dy[i]*0.5,
              head_width=0.2, head_length=0.2, fc='blue', ec='blue', alpha=0.5)

# Annotate residue index at each point
for i, r in enumerate(residues):
    plt.text(xs[i], ys[i], str(r), fontsize=7, ha='center', va='center')

plt.axis('off')
plt.tight_layout()
plt.show()
