import matplotlib.pyplot as plt
import numpy as np

# Generate the first 300 digits of Pi (excluding the initial 3) as an example
from mpmath import mp
mp.dps = 301  # Set precision to 301 digits
pi_digits = [int(d) for d in str(mp.pi)[2:302]]

# Calculate the difference ratio for Pi digits
pi_diff_ratios = [abs(pi_digits[i+1] - pi_digits[i]) / pi_digits[i] if pi_digits[i] != 0 else 0 for i in range(len(pi_digits)-1)]

# Create the horizontal line and dot positions
dot_positions = [10, 18, 28, 36, 46]  # Example positions (extend as needed)
horizontal_line_y = 0.1  # The y-coordinate for the horizontal line and dots

# Extend the dot positions to cover the range of Pi differences
dot_positions = list(range(10, len(pi_diff_ratios), 8))  # Start at 10, spaced by 8

# Plot the Pi digit difference ratios
plt.figure(figsize=(12, 6))
plt.plot(pi_diff_ratios, label="Pi Difference Ratios", color='blue', linewidth=1.5)

# Plot the horizontal line at y=0.1
plt.axhline(y=horizontal_line_y, color='orange', linestyle='--', linewidth=1, label="Horizontal Line at 0.1")

# Add dots at specified positions
for pos in dot_positions:
    if pos < len(pi_diff_ratios):  # Ensure the dot positions are within range
        plt.scatter(pos, horizontal_line_y, color='red', s=50, label="Anchor Points" if pos == 10 else "")  # Label only once

# Customize the plot
plt.title("Comparison of Pi Digit Difference Ratios with Anchor Dots")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()
plt.grid()
plt.show()
