from IPython.display import HTML
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# Triangle 1 – π-Ray collapse (1,4,1)
A1 = np.array([0, 0])
B1 = np.array([4, 0])
C1 = np.array([1, 0])

# Triangle 2 – 2,3,4 (obtuse)
A2 = np.array([0, 0])
B2 = np.array([4, 0])
C2 = np.array([2.625, 1.45237])

# Triangle 3 – 2,2,3 (isosceles)
A3 = np.array([0, 0])
B3 = np.array([3, 0])
C3 = np.array([1.5, 1.32288])

# Combine into three-phase animation
frames_per_phase = 60
total_frames = frames_per_phase * 2

# Interpolation for each triangle transition
points_A = np.concatenate([
    np.linspace(A1, A2, frames_per_phase),
    np.linspace(A2, A3, frames_per_phase)
])

points_B = np.concatenate([
    np.linspace(B1, B2, frames_per_phase),
    np.linspace(B2, B3, frames_per_phase)
])

points_C = np.concatenate([
    np.linspace(C1, C2, frames_per_phase),
    np.linspace(C2, C3, frames_per_phase)
])

# Set up figure
fig, ax = plt.subplots(figsize=(8, 5))
line, = ax.plot([], [], 'o-', lw=2, color='navy')
ax.set_xlim(-1, 5.5)
ax.set_ylim(-1, 3)
ax.set_aspect('equal')
ax.set_title('Recursive Triangle Evolution: π-Ray → 2-3-4 → 2-2-3')
ax.grid(True)

def init():
    line.set_data([], [])
    return line,

def update(frame):
    A = points_A[frame]
    B = points_B[frame]
    C = points_C[frame]
    x = [A[0], B[0], C[0], A[0]]
    y = [A[1], B[1], C[1], A[1]]
    line.set_data(x, y)
    return line,

# Animate all three
ani = animation.FuncAnimation(fig, update, frames=total_frames, init_func=init, blit=True, interval=50)

# Display in notebook
HTML(ani.to_jshtml())
