import matplotlib.pyplot as plt
import numpy as np

# The given sequence
sequence = [1, 3, 4,   3, 4, 1,   1, 1, 2,   2, 5, 5,   5, 9, 6,   6, 2, 4,   4, 6, 5,   5, 5, 4]

# Define x values
x = np.arange(len(sequence))

# Define color bands (every 3 values for triangle groups)
colors = ['#FFCCCC', '#CCE5FF', '#CCFFCC', '#FFF2CC', '#E0CCFF', '#FFD9B3', '#B3FFF0', '#FFB3E6']

# Plotting the sequence with background color bands
plt.figure(figsize=(14, 5))
for i in range(0, len(sequence), 3):
    plt.axvspan(i, i+3, color=colors[(i//3) % len(colors)], alpha=1)

plt.plot(x, sequence, marker='o', color='black', linewidth=2)
plt.title("Triadic Color-Banded Visualization of Combined Byte 1 Sequence")
plt.xlabel("Index")
plt.ylabel("Value")
plt.grid(True)
plt.xticks(x)
plt.tight_layout()
plt.show()