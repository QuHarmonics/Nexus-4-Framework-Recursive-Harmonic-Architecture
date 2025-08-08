import numpy as np
import matplotlib.pyplot as plt

# Parameters
steps = 3000000              # number of time steps
angle = 0.0               # initial angle of movement
velocity = 40.0            # base speed
beta = 0.35                # internal gain factor
correction_strength = 0.35 # feedback correction from past trend
vacuum_gain = 0.01        # how much the "vacuum" curves the path

# Initialize vacuum location (symbolic asymmetry)
x_vacuum, y_vacuum = 4096, 4096

# Initialize trajectory
x = [0]
y = [0]
angles = [angle]

for t in range(1, steps):
    # Internal curvature-based deviation
    if t > 2:
        curvature = (angles[-1] - angles[-2]) - (angles[-2] - angles[-3])
    else:
        curvature = 0
    deviation = curvature * beta

    # Recursive feedback correction
    past_trend = angles[-1] - angles[-2] if t > 1 else 0

    # Harmonic pull toward the symbolic vacuum
    dx = x_vacuum - x[-1]
    dy = y_vacuum - y[-1]
    vacuum_angle = np.arctan2(dy, dx)
    angle_difference = vacuum_angle - angles[-1]
    
    # Normalize angle difference to [-pi, pi]
    angle_difference = (angle_difference + np.pi) % (2 * np.pi) - np.pi

    # Curvature from symbolic vacuum
    attraction = angle_difference * vacuum_gain

    # Combined new angle from all feedback
    new_angle = angles[-1] + deviation - correction_strength * past_trend + attraction
    angles.append(new_angle)

    # Move in the new direction
    x.append(x[-1] + velocity * np.cos(new_angle))
    y.append(y[-1] + velocity * np.sin(new_angle))

# Plot the trajectory
plt.figure(figsize=(8, 6))
plt.plot(x, y, marker='o', markersize=2, linewidth=1.5, color='blue')
plt.plot(x_vacuum, y_vacuum, 'ro', label='Vacuum Anchor')
plt.title("Recursive Walker with Emergent Vacuum Curvature")
plt.xlabel("X Position")
plt.ylabel("Y Position")
plt.grid(True)
plt.legend()
plt.axis('equal')
plt.tight_layout()
plt.show()
