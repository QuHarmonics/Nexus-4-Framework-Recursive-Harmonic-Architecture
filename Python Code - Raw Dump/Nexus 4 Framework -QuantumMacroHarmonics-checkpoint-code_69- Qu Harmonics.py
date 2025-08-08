import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Generate a spiral flowing into a cone
def generate_cone_spiral(data_length, cone_height=10):
    theta = np.linspace(0, 4 * np.pi, data_length)  # Spiral angle
    z = np.linspace(0, cone_height, data_length)   # Height of the cone
    r = z / cone_height  # Radius decreases with height
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, z

# Visualize the cone-to-disc interaction
def plot_cone_to_disc(x, y, z, title="Cone-to-Disc Compression"):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the cone
    ax.plot(x, y, z, label="Data Flow (Cone)", color='blue')
    ax.scatter(x, y, np.zeros_like(z), label="Disc Projection (Hash)", color='orange', s=5)
    
    # Labels and titles
    ax.set_title(title)
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis (Height)")
    ax.legend()
    plt.show()

# Parameters
data_length = 512
x, y, z = generate_cone_spiral(data_length)

# Visualize the cone interacting with the disc
plot_cone_to_disc(x, y, z)
