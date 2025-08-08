# Re-import required libraries after reset
import matplotlib.pyplot as plt
import numpy as np

# Element data mapped to harmonic phases
elements = [
    "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
    "Carbon", "Nitrogen", "Oxygen", "Fluorine"
]

harmonic_phases = [
    "Note", "Echo", "Loop", "Field", "Phase Lock",
    "Harmony", "Collapse", "Compression", "Signature"
]

colors = [
    "#ff0000", "#ff6600", "#ffaa00", "#ccff00", "#66ff66",
    "#00cccc", "#3399ff", "#6633cc", "#9900cc"
]

# Spiral coordinates
theta = np.linspace(0, 4 * np.pi, len(elements))
r = np.linspace(0.5, 2.5, len(elements))
x = r * np.cos(theta)
y = r * np.sin(theta)

# Plot the spiral
plt.figure(figsize=(10, 10))
for i in range(len(elements)):
    plt.plot(x[i], y[i], 'o', markersize=18, color=colors[i])
    plt.text(x[i], y[i], f"{elements[i]}\n({harmonic_phases[i]})", ha='center', va='center', color='white', fontsize=9)

plt.plot(x, y, linestyle='--', color='gray')
plt.title("Harmonic Spiral of the First 9 Elements")
plt.axis("off")
plt.grid(False)
plt.show()
