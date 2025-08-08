import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

# === Parameters for the Spiral ===
theta = np.linspace(0, 8 * np.pi, 1000)     # Spiral angle range
z = np.linspace(0, 2, 1000)                 # Vertical range (Recursive time)
r = 0.3 + 0.05 * np.sin(9 * theta)          # Modulated radius (resonance)

x = r * np.cos(theta)
y = r * np.sin(theta)

# === Method points to highlight ===
method_indices = np.linspace(0, len(theta)-1, 9, dtype=int)

# === Set up the figure ===
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# === Initialize the scene ===
def init():
    ax.set_xlim([-0.4, 0.4])
    ax.set_ylim([-0.4, 0.4])
    ax.set_zlim([0, 2])
    ax.set_title("Nexus 2 – Harmonic Recursive Framework Spiral", fontsize=14)
    ax.set_xlabel("X (Reflection)")
    ax.set_ylabel("Y (Folding)")
    ax.set_zlabel("Z (Recursive Time)")
    return fig,

# === Update function for animation ===
def update(frame):
    ax.cla()
    ax.plot(x[:frame], y[:frame], z[:frame], color='mediumblue', linewidth=2.5)
    idx_subset = method_indices[method_indices < frame]
    ax.scatter(x[idx_subset], y[idx_subset], z[idx_subset], color='gold', s=80)
    for i, idx in enumerate(idx_subset):
        ax.text(x[idx], y[idx], z[idx], f"M{i+1}", color='black', fontsize=10, ha='center', va='center')
    ax.set_xlim([-0.4, 0.4])
    ax.set_ylim([-0.4, 0.4])
    ax.set_zlim([0, 2])
    ax.set_xlabel("X (Reflection)")
    ax.set_ylabel("Y (Folding)")
    ax.set_zlabel("Z (Recursive Time)")
    return fig,

# === Animate and save ===
spiral_animation = animation.FuncAnimation(
    fig,
    update,
    frames=np.linspace(10, 1000, 90, dtype=int),
    init_func=init,
    blit=False
)

# === Save as .mp4 ===
spiral_animation.save("Nexus2_Harmonic_Spiral_Animation.mp4", writer='ffmpeg', fps=15)

print("✅ Animation saved as Nexus2_Harmonic_Spiral_Animation.mp4")
