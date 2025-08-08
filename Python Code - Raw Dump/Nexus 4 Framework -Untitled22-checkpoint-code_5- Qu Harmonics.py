import numpy as np
import matplotlib.pyplot as plt
from hashlib import sha256
import random

# ------------------------------------------------------------------
# 1)  Helper: build unfolded numeric field from any hash-like string
# ------------------------------------------------------------------
def build_field(data: str, width: int = 128) -> np.ndarray:
    """
    Lay the (string) data row-wise into a 2-D numeric matrix.
    Non-digit characters are mapped to 0-9 by their hex value (if any),
    otherwise to 0.
    """
    # translate each character into 0-9
    digits = []
    for ch in data:
        if ch.isdigit():
            digits.append(int(ch))
        elif ch.lower() in "abcdef":
            digits.append(int(ch, 16) % 10)
        else:
            digits.append(0)

    # pad so length is multiple of width
    if len(digits) % width != 0:
        pad = width - len(digits) % width
        digits += [0] * pad

    height = len(digits) // width
    return np.array(digits, dtype=int).reshape(height, width)

# ------------------------------------------------------------------
# 2)  Simple Plinko-style probe simulator
# ------------------------------------------------------------------
def run_probes(field: np.ndarray, n_probes: int = 40, max_steps: int = 250):
    """
    Drop probes from random positions along the top edge.
    At each step a probe moves to the lowest-value neighbour among
    straight-down, down-left, down-right (simulating 'gravity').
    """
    h, w = field.shape
    paths = []

    for _ in range(n_probes):
        x = random.randrange(w)
        y = 0
        path = [(y, x)]

        for _ in range(max_steps):
            if y == h - 1:
                break  # reached bottom

            # candidate moves
            moves = []
            for dx in (-1, 0, 1):
                nx = x + dx
                ny = y + 1
                if 0 <= nx < w:
                    moves.append((field[ny, nx], ny, nx))

            # pick the move with minimal field value (strongest "gravity")
            val, y, x = min(moves, key=lambda t: t[0])
            path.append((y, x))

        paths.append(path)

    return paths

# ------------------------------------------------------------------
# 3)  Visualiser: field heatmap + drift vectors + probe paths
# ------------------------------------------------------------------
def plot_field_with_probes(field: np.ndarray, paths, show_vectors=True):
    h, w = field.shape
    plt.figure(figsize=(12, 6))
    plt.imshow(field, cmap='viridis', interpolation='nearest')
    plt.colorbar(label='Field value')

    # optional drift vectors (simple gradient)
    if show_vectors:
        # down-sampling for clarity
        step = max(1, w // 32)
        gy, gx = np.gradient(field.astype(float))
        plt.quiver(
            np.arange(0, w, step),
            np.arange(0, h, step),
            -gx[::step, ::step],
            -gy[::step, ::step],
            color='white',
            alpha=0.5,
            scale=200
        )

    # probe paths
    for path in paths:
        ys, xs = zip(*path)
        plt.plot(xs, ys, color='red', linewidth=0.8, alpha=0.7)

    plt.title("Unfolded Hash Field with Plinko Probes & Drift Vectors")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()

# ------------------------------------------------------------------
# ===  Demo  =======================================================
# ------------------------------------------------------------------
# Feel free to replace this with any 64-hex SHA-256 string or longer
sample_hash = sha256(b"Hello, harmonic breathing!").hexdigest()

# Build a toy 'reversed' / expanded string by repeated base conversion
expanded = ''.join(f"{int(ch, 16):02d}" for ch in sample_hash) * 4  # simple expansion

field = build_field(expanded, width=128)
paths = run_probes(field, n_probes=50, max_steps=field.shape[0])

plot_field_with_probes(field, paths)
