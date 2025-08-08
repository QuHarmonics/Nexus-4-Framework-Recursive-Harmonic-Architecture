#!/usr/bin/env python3
"""
spiral_vs_hexicon_fullscreen.py

Maximized side-by-side of analog spiral and hex-grid embeddings,
with no margins, huge fonts, and full use of canvas.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

# --- YOUR REAL Δ-PHASE DATA (replace with real values) ---
np.random.seed(2)
real_deltas = np.random.randn(64)
norm = (real_deltas - real_deltas.min()) / (real_deltas.max() - real_deltas.min())

# --- Analog spiral coords ---
theta = np.linspace(0, 4*np.pi, 64, endpoint=False)
radius = np.linspace(1, 5, 64)
spiral_x = radius * np.cos(theta)
spiral_y = radius * np.sin(theta)

# --- Hex-grid coords ---
def hex_spiral(n):
    coords = [(0,0)]
    layer = 1
    while len(coords) < n:
        q, r = -layer, 0
        for dq, dr in [(1,-1),(1,0),(0,1),(-1,1),(-1,0),(0,-1)]:
            for _ in range(layer):
                if len(coords) >= n:
                    return coords
                q, r = q+dq, r+dr
                coords.append((q, r))
        layer += 1
    return coords[:n]

axial = hex_spiral(64)
hex_x = np.array([np.sqrt(3)*(q + r/2) for q,r in axial])
hex_y = np.array([1.5*r for q,r in axial])

# --- Plot ---
fig = plt.figure(figsize=(24,12), dpi=200)

# Full-bleed axes for spiral (left half)
ax1 = fig.add_axes([0.01, 0.05, 0.48, 0.9])
sc1 = ax1.scatter(spiral_x, spiral_y, c=norm, cmap='viridis', s=800, edgecolors='black', linewidths=1.0)
for i,(x,y) in enumerate(zip(spiral_x,spiral_y)):
    ax1.text(x, y, f"{i:02X}", ha='center', va='center', color='white', fontsize=18, weight='bold')
ax1.set_title("Analog Spiral Embedding\nResidue 00–3F", fontsize=24, pad=20)
ax1.set_aspect('equal')
ax1.axis('off')

# Full-bleed axes for hex-grid (right half)
ax2 = fig.add_axes([0.51, 0.05, 0.48, 0.9])
R = 1.5
for i, (x, y) in enumerate(zip(hex_x, hex_y)):
    col = plt.cm.viridis(norm[i])
    hexagon = RegularPolygon(
        (x, y), numVertices=6, radius=R,
        orientation=np.pi/6, facecolor=col, edgecolor='black', linewidth=0.8
    )
    ax2.add_patch(hexagon)
    ax2.text(x, y, f"{i:02X}", ha='center', va='center', fontsize=18, color='white', weight='bold')
ax2.set_title("Digital Hexicon Grid Embedding\nResidue 00–3F", fontsize=24, pad=20)
ax2.set_aspect('equal')
ax2.axis('off')

# Big colorbar down the middle
cax = fig.add_axes([0.49, 0.25, 0.02, 0.5])
cb = plt.colorbar(sc1, cax=cax)
cb.set_label("Normalized Δ-phase", rotation=270, labelpad=25, fontsize=20)
cb.ax.tick_params(labelsize=16)

plt.show()
