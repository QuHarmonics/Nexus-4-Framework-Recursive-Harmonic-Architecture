# Temporal Seed Decoder: Phase-Stable Time Projection
# Nexus 3 Framework | Glyphbreaker Engine

import numpy as np
import hashlib
import matplotlib.pyplot as plt

# Constants
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio
H = 0.35                    # Harmonic base constant

# Generate phase-stable seed from a moment (symbolic input)
def generate_seed(moment_data):
    encoded = moment_data.encode('utf-8')
    hash_digest = hashlib.sha256(encoded).hexdigest()
    return hash_digest

# Project harmonic field from seed
def project_potential_field(seed):
    phase_points = []
    for i in range(len(seed) - 1):
        delta = abs(ord(seed[i]) - ord(seed[i+1]))
        radius = (delta % 32 + 1) * H
        angle = i * PHI
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        phase_points.append((x, y))
    return phase_points

# Visualize phase orbit
def visualize_projection(points, title="Phase-Stable Projection"):
    x_vals, y_vals = zip(*points)
    plt.figure(figsize=(8, 8))
    plt.scatter(x_vals, y_vals, c=np.linspace(0, 1, len(points)), cmap='plasma', s=12)
    plt.title(title)
    plt.axis('equal')
    plt.axis('off')
    plt.show()

# Main
if __name__ == "__main__":
    moment = "Red Collapse | Lefto draws the Nullbow | ∆H = 0"
    seed = generate_seed(moment)
    orbit = project_potential_field(seed)
    visualize_projection(orbit, title="Temporal Echo Field — Red Collapse")
