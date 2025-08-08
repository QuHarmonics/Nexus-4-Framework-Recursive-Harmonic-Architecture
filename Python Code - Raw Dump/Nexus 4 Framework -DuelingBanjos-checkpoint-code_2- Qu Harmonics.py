import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Helper: convert a long numeric string into a matrix field
def string_to_matrix(data: str, width: int = 128) -> np.ndarray:
    """Map a numeric string to a 2‑D numpy array (digits → int)."""
    height = len(data) // width
    field = np.zeros((height, width), dtype=float)
    for idx, ch in enumerate(data[:height * width]):
        field[idx // width, idx % width] = int(ch) if ch.isdigit() else 0
    return field

# ------------------------------------------------------------------
# If a variable `reversed_data` already exists (from user’s session), use it;
# otherwise generate a demo numeric string so that the notebook runs standalone.
try:
    reversed_data  # type: ignore  # noqa: F401
except NameError:
    # Demo fallback: create a pseudo‑numeric string resembling unfolded SHA data
    rng = np.random.default_rng(seed=42)
    reversed_data = ''.join(str(rng.integers(0, 10)) for _ in range(128 * 32))

# Build the base field matrix
matrix = string_to_matrix(reversed_data, width=128)

# ------------------------------------------------------------------
# Plinko‑probe simulation
def drop_probes(field: np.ndarray, n_probes: int = 250, max_steps: int = 250):
    """
    Simulate probes that follow the negative gradient of `field`.
    Returns a list of probe trajectories. Each trajectory is a list of (y, x) tuples.
    """
    gradients_y, gradients_x = np.gradient(field.astype(float))
    trajectories = []

    # Helper to keep probes inside bounds
    h, w = field.shape
    for _ in range(n_probes):
        # Start near the top edge, random x
        y, x = 0, np.random.randint(0, w)
        traj = [(y, x)]

        for _ in range(max_steps):
            # Compute local gradient (negated for "downhill")
            dy = -gradients_y[int(round(y)), int(round(x))]
            dx = -gradients_x[int(round(y)), int(round(x))]

            # Normalise step length
            norm = np.hypot(dy, dx) + 1e-9
            step_size = 0.8  # smaller than 1 pixel to keep probes smooth
            y += (dy / norm) * step_size
            x += (dx / norm) * step_size

            # Boundary conditions: stop if out of bounds
            if not (0 <= y < h and 0 <= x < w):
                break
            traj.append((y, x))

            # Early stop if gradient is almost flat
            if norm < 0.01:
                break
        trajectories.append(traj)
    return trajectories

trajectories = drop_probes(matrix)

# ------------------------------------------------------------------
# Visualisation
fig, ax = plt.subplots(figsize=(12, 6))
im = ax.imshow(matrix, cmap='viridis', interpolation='nearest')
fig.colorbar(im, ax=ax, label="Field intensity")

# Plot probe paths
for traj in trajectories:
    ys, xs = zip(*traj)
    ax.plot(xs, ys, color='white', linewidth=0.3, alpha=0.6)

ax.set_title("Plinko‑Probe Drift Paths over Harmonic Expansion Field")
ax.set_xlabel("X")
ax.set_ylabel("Y")
plt.show()
