import numpy as np
import matplotlib.pyplot as plt

# --- Minimal Plinko Universe Demo --------------------------------------------
np.random.seed(42)                  # deterministic demo

# 1. Mock "tension field" (this would normally come from SHA expansion)
field_size = 64
tension = np.random.rand(field_size, field_size)

# 2. Simple gradient‑descent plinko probe
start_x = field_size // 2           # drop from middle top
path = [(0, start_x)]
x = start_x
for y in range(field_size - 1):
    # pick next x among current, left, right based on lowest tension below
    candidates = [x]
    if x > 0:
        candidates.append(x - 1)
    if x < field_size - 1:
        candidates.append(x + 1)
    # choose the candidate with lowest tension at next row
    next_x = min(candidates, key=lambda c: tension[y + 1, c])
    path.append((y + 1, next_x))
    x = next_x

# 3. Plot the field + path
plt.figure(figsize=(6, 6))
plt.imshow(tension, cmap="viridis", origin="upper")
ys, xs = zip(*path)
plt.plot(xs, ys, color="white", lw=2, label="Plinko path")
plt.scatter([start_x], [0], color="red", s=50, label="Start")
plt.scatter([xs[-1]], [ys[-1]], color="yellow", s=50, label="End")
plt.title("Minimal Plinko Probe in Mock Tension Field")
plt.legend(loc="lower right")
plt.axis("off")
plt.show()
