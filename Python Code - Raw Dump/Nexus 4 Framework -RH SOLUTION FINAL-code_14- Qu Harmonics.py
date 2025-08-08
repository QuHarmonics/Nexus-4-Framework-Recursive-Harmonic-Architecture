import numpy as np
import matplotlib.pyplot as plt
from math import cos, pi

# Initialize Pi digits and Zeta zeros (simplified for demonstration)
pi_digits = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]
zeta_real_parts = [14.134725, 21.02204, 25.010857, 30.424876, 32.935062, 37.586178, 40.918719]

# Function to calculate gaps and ratios
def calculate_ratios(values):
    differences = np.diff(values)
    ratios = np.divide(
        differences[1:].astype(float),
        differences[:-1].astype(float),
        out=np.zeros_like(differences[1:], dtype=float),
        where=differences[:-1] != 0
    )
    return differences, ratios

# Generate triangles dynamically based on gaps and ratios
def generate_triangles(ratios, oscillations):
    triangles = []
    for i, ratio in enumerate(ratios):
        a = abs(ratio)
        b = abs(oscillations[i])
        c = np.sqrt(a**2 + b**2)  # Hypotenuse
        triangles.append((a, b, c))
    return triangles

# Plot triangles and relationships
def plot_triangles(triangles, title):
    plt.figure(figsize=(10, 8))
    for i, (a, b, c) in enumerate(triangles):
        plt.plot([0, a], [0, 0], label=f"Base {i + 1}" if i == 0 else None)  # Base
        plt.plot([0, 0], [0, b], label=f"Height {i + 1}" if i == 0 else None)  # Height
        plt.plot([a, 0], [0, b], label=f"Hypotenuse {i + 1}" if i == 0 else None)  # Hypotenuse
    plt.title(title)
    plt.xlabel("Ratio of Differences")
    plt.ylabel("Cos Oscillation")
    plt.legend()
    plt.grid(True)
    plt.show()

# Calculate cosine oscillations
def cosine_oscillations(length, frequency=pi / 4):
    return [cos(i * frequency) for i in range(length)]

# Main function to calculate and visualize
def main():
    # Calculate differences and ratios for Pi digits
    pi_differences, pi_ratios = calculate_ratios(pi_digits)
    
    # Generate cosine oscillations
    oscillations = cosine_oscillations(len(pi_ratios))
    
    # Generate triangles from ratios and oscillations
    triangles = generate_triangles(pi_ratios, oscillations)
    
    # Plot the triangles
    plot_triangles(triangles, "Dynamic Triangles from Pi Ratios and Cosine Oscillations")

    # Visualize ratios and oscillations
    plt.figure(figsize=(12, 6))
    plt.plot(pi_ratios, label="Pi Ratios", marker="o", linestyle="-", color="blue")
    plt.plot(oscillations, label="Cosine Oscillations", marker="x", linestyle="--", color="orange")
    plt.title("Ratios of Differences and Cosine Oscillations")
    plt.xlabel("Index")
    plt.ylabel("Values")
    plt.legend()
    plt.grid(True)
    plt.show()

# Run the main function
main()
