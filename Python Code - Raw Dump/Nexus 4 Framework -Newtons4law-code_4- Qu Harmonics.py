
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Sample triangle sets: (a, b, c)
triangle_sets = [(3, 4, 5), (5, 5, 6), (6, 2, 5), (2, 6, 5), (5, 2, 6)]

def triangle_coords(a, b, c):
    A = np.array([0, 0])
    B = np.array([c, 0])
    try:
        cos_angle = (a**2 + b**2 - c**2) / (2 * a * b)
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        angle_C = np.arccos(cos_angle)
        C_x = b * np.cos(angle_C)
        C_y = b * np.sin(angle_C)
        C = np.array([C_x, C_y])
    except:
        C = np.array([0, 0])
    return A, B, C

fig, ax = plt.subplots()
line, = ax.plot([], [], 'bo-', lw=2)
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_title("Matplotlib Triangle Animator")

def init():
    line.set_data([], [])
    return line,

def update(frame):
    a, b, c = triangle_sets[frame % len(triangle_sets)]
    A, B, C = triangle_coords(a, b, c)
    x = [A[0], B[0], C[0], A[0]]
    y = [A[1], B[1], C[1], A[1]]
    line.set_data(x, y)
    return line,

ani = animation.FuncAnimation(fig, update, frames=len(triangle_sets),
                              init_func=init, blit=True, repeat=True)

plt.show()
