# ─── 2-A.  parameters you can tweak easily ───────────────────────────────────
N_PROBES   = 250         # how many particles to drop
N_STEPS    = 250         # max steps per probe
STEP       = 1.0         # pixel movement per step
DROP_EDGE  = "top"       # "top", "bottom", "left", or "right"

h, w = field.shape
paths = []               # will hold a list of (x[], y[]) for each probe

# ─── 2-B.  helper to pick random spawn coord on chosen edge ──────────────────
def random_spawn(edge):
    if edge == "top":    return np.random.randint(0, w), 0
    if edge == "bottom": return np.random.randint(0, w), h - 1
    if edge == "left":   return 0, np.random.randint(0, h)
    if edge == "right":  return w - 1, np.random.randint(0, h)
    raise ValueError("edge must be top / bottom / left / right")

# ─── 2-C.  main probe loop ───────────────────────────────────────────────────
for _ in range(N_PROBES):
    x, y  = random_spawn(DROP_EDGE)
    xs, ys = [x], [y]

    for _ in range(N_STEPS):
        # sample gradient; if flat region, break
        gx, gy = GX[int(y) % h, int(x) % w], GY[int(y) % h, int(x) % w]
        gnorm  = np.hypot(gx, gy)
        if gnorm < 1e-6:      # almost zero tension → particle rests
            break
        # move downhill (steepest descent)
        x += (gx / gnorm) * STEP
        y += (gy / gnorm) * STEP
        xs.append(x); ys.append(y)

        # stop if wandered outside
        if not (0 <= x < w and 0 <= y < h):
            break
    paths.append((xs, ys))
plt.figure(figsize=(12, 6))
plt.imshow(field, cmap="viridis", origin="upper", interpolation="nearest")
plt.colorbar(label="Normalised Field Intensity")

# overlay drift (sample every N pixels for clarity)
skip = 4
plt.quiver(np.arange(0, w, skip), 
           np.arange(0, h, skip),
           GX[::skip, ::skip], 
           GY[::skip, ::skip],
           color="white",  linewidth=0.4, alpha=0.6)

# overlay Plinko probe paths
for xs, ys in paths:
    plt.plot(xs, ys, lw=0.7, alpha=0.6)

plt.title("SHA-256 Harmonic Expansion | Tension Field, Drift & Plinko Echoes")
plt.axis('off')
plt.show()
