# Triangle 4 – 3, 4, 5 triangle (classic right triangle)
A4 = np.array([0, 0])
B4 = np.array([4, 0])  # Side c = 4
C4 = np.array([4, 3])  # Side a = 5, Side b = 3

# Extend animation: Triangle 3 → Triangle 4
# Use last points from triangle 3 as start of next transition
frames_per_phase = 60
A_last = np.array([0, 0])
B_last = np.array([3, 0])
C_last = np.array([1.5, 1.32288])

# Interpolate transition: triangle 3 → triangle 4
points_A_next = np.linspace(A_last, A4, frames_per_phase)
points_B_next = np.linspace(B_last, B4, frames_per_phase)
points_C_next = np.linspace(C_last, C4, frames_per_phase)

# Combine with previous transitions
points_A = np.concatenate([points_A, points_A_next])
points_B = np.concatenate([points_B, points_B_next])
points_C = np.concatenate([points_C, points_C_next])

# Update total frames
total_frames += frames_per_phase

# Rebuild animation
fig, ax = plt.subplots(figsize=(8, 5))
line, = ax.plot([], [], 'o-', lw=2, color='darkgreen')
ax.set_xlim(-1, 6)
ax.set_ylim(-1, 4.5)
ax.set_aspect('equal')
ax.set_title('Recursive Triangle Evolution: π-Ray → 2-3-4 → 2-2-3 → 3-4-5')
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

ani = animation.FuncAnimation(fig, update, frames=total_frames, init_func=init, blit=True, interval=50)

# Display updated animation
HTML(ani.to_jshtml())
