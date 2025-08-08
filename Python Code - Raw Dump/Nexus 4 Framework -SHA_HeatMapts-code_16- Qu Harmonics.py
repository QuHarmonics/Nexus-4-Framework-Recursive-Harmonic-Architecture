import hashlib
import numpy as np
import matplotlib.pyplot as plt

# ——— PARAMETERS ———
initial_state = 0
iterations    = 300
MOD_RESIDUE   = 256   # or 64 if using mod-64 residues

# ——— SHA‐256 / Δ‐PHASE WALK ———
states = [initial_state]
for _ in range(iterations):
    s = states[-1]
    # 1) Compute SHA‐256 of the single‐byte state
    h = hashlib.sha256(bytes([s % 256])).digest()
    # 2) Convert to integer, normalize
    h_int = int.from_bytes(h, 'big')
    # 3) Projection b = H(s) mod MOD_RESIDUE
    b = h_int % MOD_RESIDUE
    # 4) Projection a = s           mod MOD_RESIDUE
    a = s % MOD_RESIDUE
    # 5) Δ
    d = (b - a) % MOD_RESIDUE
    # 6) Next state
    s_next = (s + d) % MOD_RESIDUE
    states.append(s_next)

# Trim to length T
T = len(states) - 1
states      = np.array(states)
next_states = states[1:]
states      = states[:-1]

# ——— BUILD DATA ARRAYS ———
# delta_residue[t] = how much we jumped (0…MOD_RESIDUE−1)
delta_residue = (next_states - states) % MOD_RESIDUE

# Recompute hashes & h_norm
h_norm = []
for s in states:
    h = hashlib.sha256(bytes([s % 256])).digest()
    h_int = int.from_bytes(h, 'big')
    h_norm.append(h_int / 2**256)
h_norm = np.array(h_norm)

# Build bit array of the residue (width = ceil(log2(MOD_RESIDUE)))
# Example for mod-256 we use 8 bits (7…0)
width = int(np.ceil(np.log2(MOD_RESIDUE)))
bit_array = ((delta_residue[:, None] & (1 << np.arange(width-1, -1, -1))) > 0).astype(int)

# msb_flip[t] = did the top bit flip between t and t+1?
# we need one fewer entry than states→ use the change in bit_array
msb_flip = np.abs(np.diff(bit_array[:, 0]))  # bit_array[:,0] is MSB
# pad to length T
msb_flip = np.concatenate([msb_flip, [0]])   # no flip at final dummy step

# ——— PLOT 1: Scatter Δ vs MSB Flip colored by h_norm ———
plt.figure(figsize=(8,4))
plt.scatter(delta_residue, msb_flip, c=h_norm, cmap='RdYlBu', edgecolor='k', alpha=0.7)
plt.colorbar(label='Normalized Hash $H_{norm}$')
plt.xlabel('Residue Δ')
plt.ylabel('MSB Flip (0=no, 1=yes)')
plt.title('Residue Δ vs. MSB Flip — Colored by $H_{norm}$')
plt.yticks([0,1])
plt.tight_layout()
plt.show()

# ——— PLOT 2: Flip‐rate histogram per residue‐bin ———
bins = np.arange(0, MOD_RESIDUE+1) - 0.5
flip_count, _ = np.histogram(delta_residue[msb_flip==1], bins=bins)
total_count, _ = np.histogram(delta_residue,        bins=bins)
rate = flip_count / np.maximum(total_count, 1)

centers = (bins[:-1] + bins[1:]) / 2
plt.figure(figsize=(8,3))
plt.bar(centers, rate, width=1.0, edgecolor='k')
plt.xlabel('Residue Δ')
plt.ylabel('P(MSB Flip)')
plt.title('Flip Probability per Residue Δ')
plt.tight_layout()
plt.show()
