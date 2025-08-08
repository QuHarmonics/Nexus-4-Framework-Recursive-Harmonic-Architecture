import numpy as np
import matplotlib.pyplot as plt
import random, math

# ------------------------------------------------------------
# 1.  HELPER – build a numeric matrix from an input string
# ------------------------------------------------------------
def string_to_field(s: str, width: int = 128) -> np.ndarray:
    """
    Map an arbitrary alphanumeric string to a 2-D numeric field.

    - `width` sets the number of columns.
    - Characters 0-9 are mapped to their int value (0-9).
    - Hex letters a-f / A-F are mapped 10-15.
    - Anything else is mapped to 0.
    """
    # translate char → int
    def c2i(c):
        if c.isdigit():
            return int(c)
        if c.lower() in "abcdef":
            return 10 + "abcdef".index(c.lower())
        return 0
    
    vals = [c2i(c) for c in s]
    
    # pad so len is multiple of width
    remainder = len(vals) % width
    if remainder:
        vals.extend([0]*(width - remainder))
        
    field = np.array(vals, dtype=float).reshape(-1, width)
    return field

# ------------------------------------------------------------
# 2.  GRADIENT (“tension”) of the field
# ------------------------------------------------------------
def tension_field(field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Return x-gradient, y-gradient (finite differences).
    """
    # forward differences, pad last row/col with 0
    dx = np.zeros_like(field)
    dy = np.zeros_like(field)
    dx[:, :-1] = field[:, 1:] - field[:, :-1]
    dy[:-1, :] = field[1:, :] - field[:-1, :]
    return dx, dy

# ------------------------------------------------------------
# 3.  PLINKO PROBE SIMULATION
# ------------------------------------------------------------
def drop_probe(field: np.ndarray,
               start_col: int | None = None,
               max_steps: int = 500) -> list[tuple[int, int]]:
    """
    Simulate a probe moving through the tension field.
    
    The probe starts at the *top* row (y = 0). At each step it moves
    1 pixel in the direction of steepest descent in |gradient|.
    Stops when it leaves the field or after max_steps.
    """
    rows, cols = field.shape
    if start_col is None:
        start_col = random.randint(0, cols-1)
        
    x, y = start_col, 0
    path = [(y, x)]
    
    dx, dy = tension_field(field)
    
    for _ in range(max_steps):
        # if probe left the field → stop
        if not (0 <= x < cols and 0 <= y < rows):
            break
        
        # compute local gradient
        gx = dx[y, x]
        gy = dy[y, x]
        mag = math.hypot(gx, gy)
        
        # if no gradient → random tiny drift downward
        if mag == 0:
            y += 1
        else:
            # normalize, move one pixel toward +|gradient| direction
            # note dy is “down” direction already.
            step_x = int(round(gx / mag))
            step_y = int(round(gy / mag))
            # ensure at least moves downward a bit so we eventually exit
            if step_y == 0:
                step_y = 1 if gy >= 0 else -1
            x += step_x
            y += step_y
        
        path.append((y, x))
        if y >= rows or y < 0:
            break
    
    return path

# ------------------------------------------------------------
# 4.  VISUALISE FIELD + PROBES
# ------------------------------------------------------------
def plot_field_with_probes(field: np.ndarray,
                           probes: list[list[tuple[int,int]]],
                           title: str = "Plinko Probe Paths"):
    plt.figure(figsize=(12,6))
    plt.imshow(field, cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Field value")
    
    colors = plt.cm.rainbow(np.linspace(0,1,len(probes)))
    for path, color in zip(probes, colors):
        ys, xs = zip(*path)
        plt.plot(xs, ys, color=color, linewidth=2, alpha=0.8)
    
    plt.title(title)
    plt.gca().invert_yaxis()   # so y=0 (top) at top of plot
    plt.show()

# ------------------------------------------------------------
# 5.  DEMO — use a SHA-256 like hex string as our field
# ------------------------------------------------------------
sha_hex = "5e36e1b0204be9a4e3d4ed56b05058db73a1bc639027959a34ecb8dcb7ec5c91"
field = string_to_field(sha_hex*4, width=128)   # repeat hex to fill some rows

# generate a few probes
probe_paths = [drop_probe(field, start_col=random.randint(0,127)) for _ in range(8)]

plot_field_with_probes(field, probe_paths,
                       title="SHA-256 Harmonic Expansion – Plinko Probe Visualization")

