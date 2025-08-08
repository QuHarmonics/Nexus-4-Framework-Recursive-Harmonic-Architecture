#!/usr/bin/env python3
"""
prime_attractor_histogram.py

Histogram of MSB‐flip nearest SHA-256 constants in the index range 34–46,
colored by their arithmetic progression class mod 9 (prime-rich vs. prime-poor).
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

# --- Replace with your actual nearest_js list ---
# nearest_js: list of K_j indices (0–63) nearest to residue at each MSB flip
# For demonstration, we'll simulate some data:
np.random.seed(1)
nearest_js = list(np.random.choice(range(64), size=200, p=None))

# Filter to j in [34..46]
window = list(range(34, 47))
filtered = [j for j in nearest_js if j in window]
counts = Counter(filtered)

# Prepare data for plotting
js = sorted(window)
heights = [counts[j] for j in js]

# Classify each j by j mod 9
# Prime-poor classes: {0,3,6}, prime-rich: {1,2,4,5,7,8}
prime_poor = {0,3,6}
bar_colors = []
for j in js:
    cls = j % 9
    if cls in prime_poor:
        bar_colors.append('lightgray')
    else:
        bar_colors.append('steelblue')

# Plot
plt.figure(figsize=(8,4))
plt.bar(js, heights, color=bar_colors, edgecolor='black')
plt.xticks(js)
plt.xlabel('SHA-256 Constant Index j')
plt.ylabel('MSB-Flip Events Count')
plt.title('MSB-Flip Nearest Constants in j=34…46\n(blue=prime-rich mod 9, gray=prime-poor mod 9)')
plt.tight_layout()
plt.show()
