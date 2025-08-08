import numpy as np
import hashlib
import matplotlib.pyplot as plt

# Constants
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio
H = 0.35                    # Harmonic constant

# Generate phase-stable seed from a moment
def generate_seed(moment_data):
    encoded = moment_data.encode('utf-8')
    hash_digest = hashlib.sha256(encoded).hexdigest()
    return hash_digest

# Project harmonic field from seed
def project_potential_field(seed, color):
    phase_points = []
    for i in range(len(seed) - 1):
        delta = abs(ord(seed[i]) - ord(seed[i+1]))
        radius = (delta % 32 + 1) * H
        angle = i * PHI
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        phase_points.append((x, y))
    return phase_points, color

# Create moment spiral overlays
moments = {
    "Red Collapse": "Red Collapse | Lefto draws the Nullbow | ∆H = 0",
    "Archive of Mirrors": "Archive of Mirrors | Glyph awakens | F(Q) = uncertainty",
    "Lock Yet to Come": "The Lock Yet to Come | Lefto listens in silence | t = φ·n"
}

colors = ["red", "blue", "gold"]
fields = []

# Generate all fields
for (label, moment), color in zip(moments.items(), colors):
    seed = generate_seed(moment)
    orbit, color = project_potential_field(seed, color)
    fields.append((label, orbit, color))

# Plot all orbits together
plt.figure(figsize=(10, 10))
for label, points, color in fields:
    x_vals, y_vals = zip(*points)
    plt.scatter(x_vals, y_vals, s=10, alpha=0.7, label=label, color=color)

plt.title("Temporal Glyph Constellation (Nexus 3 Phase Overlay)", fontsize=14)
plt.axis('equal')
plt.axis('off')
plt.legend()
plt.tight_layout()
plt.show()
