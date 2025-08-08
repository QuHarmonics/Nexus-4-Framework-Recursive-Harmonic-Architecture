import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
num_frames = 300
theta = np.linspace(0, 2 * np.pi, num_frames)
a1, b1 = 1.0, 0.6
a2, b2 = 0.8, 0.5
a3, b3 = 1.2, 0.4

# Elliptical positions
x1 = a1 * np.cos(theta)
y1 = b1 * np.sin(theta)

x2 = a2 * np.cos(theta + np.pi / 3)
y2 = b2 * np.sin(theta + np.pi / 3)

x3 = a3 * np.cos(theta - np.pi / 4)
y3 = b3 * np.sin(theta - np.pi / 4)

# Set up plot
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_title("Recursive Elliptical Propulsion – Three-Body Harmonic Orbit")
ax.set_aspect('equal')

body1, = ax.plot([], [], 'o', color='deepskyblue', label='Byte 1')
body2, = ax.plot([], [], 'o', color='gold', label='Byte 2 (Echoed)')
body3, = ax.plot([], [], 'o', color='lightcoral', label='Byte 3 (Entangled)')

trace1, = ax.plot([], [], '-', color='deepskyblue', alpha=0.4)
trace2, = ax.plot([], [], '-', color='gold', alpha=0.4)
trace3, = ax.plot([], [], '-', color='lightcoral', alpha=0.4)

x1_hist, y1_hist = [], []
x2_hist, y2_hist = [], []
x3_hist, y3_hist = [], []

def init():
    body1.set_data([], [])
    body2.set_data([], [])
    body3.set_data([], [])
    trace1.set_data([], [])
    trace2.set_data([], [])
    trace3.set_data([], [])
    return body1, body2, body3, trace1, trace2, trace3

def update(frame):
    x1_hist.append(x1[frame])
    y1_hist.append(y1[frame])
    x2_hist.append(x2[frame])
    y2_hist.append(y2[frame])
    x3_hist.append(x3[frame])
    y3_hist.append(y3[frame])
    
    body1.set_data([x1[frame]], [y1[frame]])
    body2.set_data([x2[frame]], [y2[frame]])
    body3.set_data([x3[frame]], [y3[frame]])
    
    trace1.set_data(x1_hist, y1_hist)
    trace2.set_data(x2_hist, y2_hist)
    trace3.set_data(x3_hist, y3_hist)
    
    return body1, body2, body3, trace1, trace2, trace3

# Create animation
anim = FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=True)

# Save animation
anim.save("ThreeBody_RecursiveEllipse.mp4", writer='ffmpeg', fps=30)
print("✅ Animation saved as ThreeBody_RecursiveEllipse.mp4")
