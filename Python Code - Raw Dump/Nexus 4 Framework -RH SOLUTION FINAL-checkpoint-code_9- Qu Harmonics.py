# Generate a harmonic waveform based on π and project it onto the circle
# We use the first few digits of π to create a harmonic oscillation
pi_digits = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]
normalized_pi_digits = [digit / 10 for digit in pi_digits]  # Normalize for scaling

# Create angles for π digits in the circular space
pi_angles = np.linspace(0, 2 * np.pi, len(normalized_pi_digits))

# Generate circular coordinates for π-derived harmonic oscillation
pi_x_coords = [np.cos(angle) * magnitude for angle, magnitude in zip(pi_angles, normalized_pi_digits)]
pi_y_coords = [np.sin(angle) * magnitude for angle, magnitude in zip(pi_angles, normalized_pi_digits)]

# Plot the Zeta zeros and the π-derived harmonic oscillation on the same circular plot
fig, ax = plt.subplots(figsize=(8, 8))

# Draw the unit circle
circle = plt.Circle((0, 0), 1, color='lightgray', fill=False, linewidth=2)
ax.add_artist(circle)

# Plot Zeta zero projections
ax.scatter(x_coords, y_coords, color='blue', label='Zeta Zeros', s=50)

# Annotate Zeta zeros
for i, (x, y) in enumerate(zip(x_coords, y_coords)):
    ax.text(x * 1.1, y * 1.1, f"Z{i+1}", color='darkblue', fontsize=9)

# Plot π harmonic oscillation
ax.plot(pi_x_coords, pi_y_coords, color='red', label='π Harmonic Oscillation', linewidth=2, linestyle='--')

# Aesthetic adjustments
ax.set_xlim([-1.2, 1.2])
ax.set_ylim([-1.2, 1.2])
ax.axhline(0, color='black', linewidth=0.5, linestyle='--')
ax.axvline(0, color='black', linewidth=0.5, linestyle='--')
ax.set_aspect('equal', adjustable='datalim')
ax.set_title("Circular Projection of Zeta Zeros with π Harmonic Oscillation", fontsize=14)
ax.legend(loc='upper right')

plt.show()
