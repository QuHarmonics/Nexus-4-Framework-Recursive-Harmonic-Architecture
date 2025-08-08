from mpl_toolkits.mplot3d import Axes3D

# 3D Visualization of Lattice Growth
def visualize_3d_lattice(base, iterations=500, initial_value=1):
    x, y, z = [], [], []
    value = initial_value
    for n in range(iterations):
        x.append(value * np.cos(n * np.pi / iterations))
        y.append(value * np.sin(n * np.pi / iterations))
        z.append(value)
        value *= (-base) / (n + 1)
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x, y, z, label=f"Base {base} Lattice")
    ax.set_title("3D Lattice Growth")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.legend()
    plt.show()

# Visualize for Base -2
visualize_3d_lattice(-2)
