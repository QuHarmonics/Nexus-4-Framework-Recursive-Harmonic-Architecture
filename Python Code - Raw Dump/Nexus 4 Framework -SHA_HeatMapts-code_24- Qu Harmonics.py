#!/usr/bin/env python3
"""
hexicon_grid_improved.py

Larger, centered hexagonal grid of 6-bit residues with color‐coded Δ-phase.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
import hashlib

# --- STEP 1: Generate hex coordinates for 64 residues in a roughly hexagonal layout ---
# Use “spiral” coordinates so they form a nice compact hexagon
def hex_spiral(n):
    """Return first n axial hex coords in a spiral around (0,0)."""
    coords = [(0, 0)]
    layer = 1
    while len(coords) < n:
        # start at (−layer, 0)
        q, r = -layer, 0
        # six directions
        directions = [(1, -1), (1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1)]
        for direction in directions:
            for _ in range(layer):
                if len(coords) >= n: 
                    return coords
                q, r = q + direction[0], r + direction[1]
                coords.append((q, r))
        layer += 1
    return coords[:n]

coords = hex_spiral(64)

# Convert axial to Cartesian (pointy‐top)
coords_xy = []
for q, r in coords:
    x = np.sqrt(3) * (q + r/2)
    y = 1.5 * r
    coords_xy.append((x, y))
coords_xy = np.array(coords_xy)

# --- STEP 2: Compute toy Δ-phase for each residue (replace with your real data) ---
K0_mod64 = int(hashlib.sha256(b'').digest()[0]) & 0x3F
deltas = np.array([(((K0_mod64 - r) % 64) - 32) for r in range(64)])
norm = (deltas - deltas.min()) / (deltas.max() - deltas.min())

# --- STEP 3: Plot with larger figure, radius, and centered axes ---
fig, ax = plt.subplots(figsize=(10,10))

# hexagon radius
R = 1.0

for i, (x, y) in enumerate(coords_xy):
    color = plt.cm.viridis(norm[i])
    hexagon = RegularPolygon(
        (x, y), numVertices=6, radius=R, 
        orientation=np.pi/6, facecolor=color, edgecolor='gray', linewidth=0.5
    )
    ax.add_patch(hexagon)
    ax.text(x, y, f"{i:02X}", ha='center', va='center', fontsize=9, color='white')

# compute limits to center the grid
xs, ys = coords_xy[:,0], coords_xy[:,1]
pad = R * 1.5
ax.set_xlim(xs.min() - pad, xs.max() + pad)
ax.set_ylim(ys.min() - pad, ys.max() + pad)

ax.set_aspect('equal')
ax.axis('off')

cbar = fig.colorbar(
    plt.cm.ScalarMappable(cmap='viridis'),
    ax=ax, fraction=0.046, pad=0.04
)
cbar.set_label(r'Normalized $\Delta(r)$', rotation=270, labelpad=15)

plt.title("Hexicon Grid ( Improved Layout & Size )", fontsize=16)
plt.tight_layout()
plt.show()
