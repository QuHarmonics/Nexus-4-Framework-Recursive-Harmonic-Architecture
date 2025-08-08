# Re-import required libraries after code execution environment reset
import numpy as np
import matplotlib.pyplot as plt
from mpmath import mp

# Set high precision for pi digit extraction
mp.dps = 500  # decimal places

# Grid size and ray parameters
L = 8
num_bounces = 400

origin = np.array([0, 0])
direction = np.array([3.0, 0.35])
direction = direction / np.linalg.norm(direction)

# Compute ray path
points = [origin.copy()]
position = origin.copy()
velocity = direction.copy()

visited = set()
echo_mass_over_time = []

for _ in range(num_bounces):
    next_position = position + velocity
    for i in range(2):
        if next_position[i] < 0 or next_position[i] > L:
            velocity[i] = -velocity[i]
            next_position = position + velocity
    position = next_position
    ix, iy = int(np.floor(position[0])), int(np.floor(position[1]))
    visited.add((ix, iy))
    echo_mass_over_time.append(len(visited))
    points.append(position.copy())

# Get 400 hex digits of pi and convert to 2D coordinates
pi_hex_digits = mp.nstr(mp.pi, n=500, strip_zeros=False)
pi_digits = [int(ch, 16) for ch in pi_hex_digits if ch in '0123456789abcdef'[:16]][2:402]
pi_coords = [(pi_digits[i] % L, pi_digits[i + 1] % L) for i in range(0, len(pi_digits)-1, 2)]

# Convert visited positions to array for plotting
visited_x, visited_y = zip(*visited)
pi_x, pi_y = zip(*pi_coords)

# Plot lattice overlay of echo mass vs BBP-aligned digits
plt.figure(figsize=(8, 8))
plt.scatter(visited_x, visited_y, color='orange', label='Echo Points', alpha=0.6)
plt.scatter(pi_x, pi_y, color='blue', marker='x', label='π Digit Projections', alpha=0.5)
plt.title("Echo Grid vs π-ray Projections (Vector [3.0, 0.35])")
plt.xlabel("X")
plt.ylabel("Y")
plt.xlim([0, L])
plt.ylim([0, L])
plt.grid(True)
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')

plt.tight_layout()
plt.show()
