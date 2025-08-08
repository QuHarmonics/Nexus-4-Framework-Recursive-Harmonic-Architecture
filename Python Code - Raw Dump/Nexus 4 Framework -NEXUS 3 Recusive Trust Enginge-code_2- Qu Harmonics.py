import numpy as np
import matplotlib.pyplot as plt

def plot_lean(data, width):
    """Plots the data string as a 2D grid and shows density."""
    height = len(data) // width
    matrix = np.zeros((height, width))
    
    for i, char in enumerate(data[:height*width]):
        matrix[i//width, i%width] = int(char) if char.isdigit() else 0  # treat non-digits as 0
    
    plt.figure(figsize=(12, 6))
    plt.imshow(matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    plt.title("Harmonic Expansion Field - Drift Visualization")
    plt.show()
    
    return matrix

# Example:
matrix = plot_lean(reversed_data, width=99)
