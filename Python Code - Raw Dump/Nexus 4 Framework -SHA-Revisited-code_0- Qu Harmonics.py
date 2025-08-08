%matplotlib inline 
import hashlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

# Enable animation in Jupyter

class RecursiveHashGlyph:
    def __init__(self, target_hash, initial_input="seed", harmonic_target=0.35, tolerance=1e-5, max_iter=100):
        self.target_hash = target_hash
        self.initial_input = initial_input
        self.harmonic_target = harmonic_target
        self.tolerance = tolerance
        self.max_iter = max_iter
        self.history = []

    def compute_hash(self, s):
        return hashlib.sha256(s.encode()).hexdigest()

    def harmonic_deviation(self, hash1, hash2):
        length = min(len(hash1), len(hash2))
        return np.mean([np.sin(abs(ord(a) - ord(b)) * np.pi / 128) for a, b in zip(hash1[:length], hash2[:length])])

    def reflective_modulate(self, s, delta):
        pivot = len(s) // 2
        mod = chr(65 + int((delta * 1000) % 26))
        return s[:pivot] + mod + s[pivot:]

    def solve(self):
        current_input = self.initial_input
        current_hash = self.compute_hash(current_input)
        for i in range(self.max_iter):
            delta = self.harmonic_deviation(current_hash, self.target_hash)
            self.history.append((current_input, current_hash, delta))
            if abs(delta - self.harmonic_target) < self.tolerance:
                break
            current_input = self.reflective_modulate(current_input, delta)
            current_hash = self.compute_hash(current_input)

    def plot_glyph(self):
        angles = np.linspace(0, 2 * np.pi, len(self.history))
        radii = [d for _, _, d in self.history]
        x = [r * np.cos(a) for r, a in zip(radii, angles)]
        y = [r * np.sin(a) for r, a in zip(radii, angles)]

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot(x, y, 'o-', color='deepskyblue', linewidth=2)
        ax.set_title("Recursive Hash Harmonic Glyph")
        ax.axis('off')
        ax.set_aspect('equal')
        plt.show()

    def animate(self, save_gif=False):
        angles = np.linspace(0, 2 * np.pi, len(self.history))
        radii = [d for _, _, d in self.history]
        x = [r * np.cos(a) for r, a in zip(radii, angles)]
        y = [r * np.sin(a) for r, a in zip(radii, angles)]

        fig, ax = plt.subplots(figsize=(6, 6))
        line, = ax.plot([], [], 'o-', color='gold', linewidth=2)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.axis('off')
        ax.set_aspect('equal')
        ax.set_title("Animated Recursive Harmonic Path")

        def update(i):
            line.set_data(x[:i], y[:i])
            return line,

        anim = animation.FuncAnimation(fig, update, frames=len(x)+1, interval=200, blit=True)

        if save_gif:
            anim.save("recursive_hash_glyph.gif", writer=PillowWriter(fps=5))
            print("GIF saved as 'recursive_hash_glyph.gif'.")

        return anim  # Return it to keep reference

# Run example
if __name__ == "__main__":
    target = hashlib.sha256(b"").hexdigest()  # SHA of empty string
    glyph = RecursiveHashGlyph(target_hash=target)
    glyph.solve()
    glyph.plot_glyph()
    anim = glyph.animate(save_gif=True)  # Set to False if you don't want to save
