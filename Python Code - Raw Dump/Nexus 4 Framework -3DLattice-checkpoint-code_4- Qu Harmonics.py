def inspect_lattice_center(lattice):
    """
    Inspect the center region of the lattice and visualize its values.
    """
    lattice_size = lattice.shape[0]
    center = lattice_size // 2  # Center index

    # Extract a smaller sub-cube around the center
    sub_cube_size = 5
    center_slice = slice(center - sub_cube_size, center + sub_cube_size)
    sub_cube = lattice[center_slice, center_slice, center_slice]

    # Visualize the sub-cube
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = np.nonzero(sub_cube)
    values = sub_cube[x, y, z]
    scatter = ax.scatter(x, y, z, c=values, cmap='viridis', s=20)

    ax.set_title("Center Region of the Lattice", fontsize=14)
    ax.set_xlabel("X-axis", fontsize=10)
    ax.set_ylabel("Y-axis", fontsize=10)
    ax.set_zlabel("Z-axis", fontsize=10)
    plt.colorbar(scatter, ax=ax, label="Harmonic Values")
    plt.show()

# Call the function
inspect_lattice_center(lattice)
