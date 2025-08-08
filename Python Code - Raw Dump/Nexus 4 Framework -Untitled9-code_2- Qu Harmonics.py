import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import random

# ------------------------------------------------------------------
# 1.  INPUT  (You can swap this out with ANY unfolded / expanded hash
#             string you like – just keep it as a long 0‑9 string)
# ------------------------------------------------------------------
reversed_data = (
    "5143614111020411149104143134141356110505813117"
    "..."  # ← snipped for brevity – feel free to paste the full string!
)
# choose a width that divides the total nicely; 128 is a good default
WIDTH = 128

# ------------------------------------------------------------------
# 2.  BUILD A 2‑D “TENSION FIELD” MATRIX FROM THE STRING
# ------------------------------------------------------------------
length = len(reversed_data)
HEIGHT = length // WIDTH
field = np.zeros((HEIGHT, WIDTH), dtype=float)

for i, ch in enumerate(reversed_data[:HEIGHT * WIDTH]):
    # simple numeric value: 0‑9 → 0‑9, anything else → 0  (space for hex/letters)
    field[i // WIDTH, i % WIDTH] = float(ch) if ch.isdigit() else 0.0

# ------------------------------------------------------------------
# 3.  DRIFT (GRADIENT) VECTORS  –  ∇field  ≈   local turbulence flow
# ------------------------------------------------------------------
dx = np.zeros_like(field)
dy = np.zeros_like(field)

dx[:, :-1] = field[:, 1:] - field[:, :-1]
dy[:-1, :] = field[1:, :] - field[:-1, :]

# ------------------------------------------------------------------
# 4.  PLINKO PROBE SIMULATION
#     • Drop N probes from random positions on the top row
#     • At each step, move one pixel in the direction of the steepest
#       descent (largest negative gradient magnitude)
# ------------------------------------------------------------------
N_PROBES   = 30        # number of probes
MAX_STEPS  = 300       # steps per probe
probe_paths = []

for _ in range(N_PROBES):
    x = random.randint(0, WIDTH - 1)   # starting column
    y = 0                              # start at top row
    path = [(y, x)]

    for _ in range(MAX_STEPS):
        if y >= HEIGHT - 1:            # reached bottom
            break
        # look at the 3 downward neighbours: (y+1, x-1), (y+1, x), (y+1, x+1)
        choices = []
        for dx_step in (-1, 0, 1):
            nx = x + dx_step
            ny = y + 1
            if 0 <= nx < WIDTH and ny < HEIGHT:
                # effective “potential” is the field value plus gradient magnitude
                gmag = abs(dx[ny, nx]) + abs(dy[ny, nx])
                choices.append((field[ny, nx] + gmag, nx, ny))
        # choose the smallest “potential” descent
        choices.sort(key=lambda t: t[0])
        _, x, y = choices[0]
        path.append((y, x))
    probe_paths.append(path)

# ------------------------------------------------------------------
# 5.  VISUALIZATION
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(14, 6))

# base heatmap of the tension field
im = ax.imshow(field, cmap=cm.viridis, interpolation='nearest', aspect='auto')
fig.colorbar(im, ax=ax, label="Field Intensity")

# drift vectors (quiver)
step = 4  # subsample arrows so the plot isn't too busy
Y, X = np.mgrid[0:HEIGHT:step, 0:WIDTH:step]
ax.quiver(
    X, Y,
    dx[::step, ::step], dy[::step, ::step],
    color="white", pivot='mid', scale=100, headwidth=3, headlength=4
)

# overlay plinko probe paths
for path in probe_paths:
    yy, xx = zip(*path)
    ax.plot(xx, yy, lw=0.8, alpha=0.8)

ax.set_title("Harmonic Expansion Field  •  Drift Vectors  •  Plinko Probes")
ax.set_xlabel("X")
ax.set_ylabel("Y")
plt.tight_layout()
plt.show()
