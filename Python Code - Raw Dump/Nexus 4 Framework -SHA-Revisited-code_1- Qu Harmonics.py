import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from hashlib import sha256
from matplotlib.animation import PillowWriter
from IPython.display import HTML

# --- Coulomb-like Collapse Spiral Simulator ---

def harmonic_potential(delta_phi):
    return 1 - np.cos(delta_phi)

def generate_spiral(seed, steps=32):
    hash_value = sha256(seed.encode()).hexdigest()
    radii = []
    angles = []

    for i in range(steps):
        byte_pair = hash_value[2 * i: 2 * i + 2]
        byte_val = int(byte_pair, 16)
        delta_phi = (byte_val / 255.0) * np.pi  # Normalize to [0, pi]
        V = harmonic_potential(delta_phi)
        r = V  # Inverse collapse strength
        angle = (2 * np.pi / steps) * i
        radii.append(r)
        angles.append(angle)

    x = [r * np.cos(a) for r, a in zip(radii, angles)]
    y = [r * np.sin(a) for r, a in zip(radii, angles)]
    return x, y, radii

def animate_spiral(x, y, radii):
    fig, ax = plt.subplots(figsize=(6, 6))
    line, = ax.plot([], [], 'o-', color='crimson', linewidth=2)
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.set_aspect('equal')
    ax.set_title("Coulomb-Inspired SHA Collapse Spiral")
    ax.axis('off')

    def update(i):
        line.set_data(x[:i], y[:i])
        return line,

    ani = animation.FuncAnimation(fig, update, frames=len(x)+1, interval=100, blit=True)
    return ani

# Example usage:
if __name__ == "__main__":
    seed_input = "light"  # You can change this to any phrase
    x_vals, y_vals, r_vals = generate_spiral(seed_input)
    anim = animate_spiral(x_vals, y_vals, r_vals)
    HTML(anim.to_jshtml())
