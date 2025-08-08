# %matplotlib notebook # Comment this out - we will use a different display method
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import matplotlib.animation as animation # Ensure animation is imported

# Import HTML display utility from IPython
from IPython.display import HTML

# Triangle 1 (Pi-Ray degenerate triangle): All points on a straight line
A1 = np.array([0, 0])
B1 = np.array([4, 0])
C1 = np.array([1, 0]) # Linear collapse triangle (degenerate)

# Triangle 2 (Emergent triangle from 2, 3, 4 sides)
A2 = np.array([0, 0])
B2 = np.array([4, 0])
C2 = np.array([2.625, 1.45237]) # As per triangle calculator output

# Animation setup
frames = 60
points_A = np.linspace(A1, A2, frames)
points_B = np.linspace(B1, B2, frames) # B point is static as B1=B2
points_C = np.linspace(C1, C2, frames)

# Set up plot
fig, ax = plt.subplots(figsize=(8, 5))
line, = ax.plot([], [], 'o-', lw=2, color='blue')
ax.set_xlim(-1, 5.5)
ax.set_ylim(-1, 3)
ax.set_aspect('equal')
ax.set_title('Recursive Harmonic Triangle Transition\n(π-Ray → Emergent 2-3-4 Triangle)')
ax.grid(True)

# Close the figure to prevent it from being displayed as a static plot before the animation
plt.close(fig)

# Initialization function
def init():
    line.set_data([], [])
    return line,

# Update function
def update(frame):
    A = points_A[frame]
    B = points_B[frame]
    C = points_C[frame]
    x = [A[0], B[0], C[0], A[0]]
    y = [A[1], B[1], C[1], A[1]]
    line.set_data(x, y)
    return line,

# Run animation
# Need to ensure a writer like 'ffmpeg' or 'imagemagick' is available or use the html5 backend
# Use the HTML5 backend, which requires no external writers but relies on browser capabilities
ani = animation.FuncAnimation(fig, update, frames=frames, init_func=init, blit=True, interval=50)

# Show the animation using HTML5 video
# This method generates HTML for an HMTL5 video tag
# HTML(ani.to_jshtml()) # Alternative using JS, might still hit IPython if JS is the issue
HTML(ani.to_html5_video()) # Renders as HTML5 video tag

# plt.show() # This is no longer needed with the HTML display method