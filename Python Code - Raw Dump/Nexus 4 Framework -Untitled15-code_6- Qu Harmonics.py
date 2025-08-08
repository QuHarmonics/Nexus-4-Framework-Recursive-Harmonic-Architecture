import matplotlib.pyplot as plt
import numpy as np

# Generate the first 300 digits of Pi (excluding the initial 3) as an example
from mpmath import mp
mp.dps = 301  # Set precision to 301 digits
pi_digits = [int(d) for d in str(mp.pi)[2:302]]

# Calculate the difference ratio for Pi digits
pi_diff_ratios = [abs(pi_digits[i+1] - pi_digits[i]) / pi_digits[i] if pi_digits[i] != 0 else 0 for i in range(len(pi_digits)-1)]

# Create the repeating anchor pattern (10, 8, 10)
anchor_pattern = []
current_value = 0
for i in range(len(pi_diff_ratios)):
    if i % 3 == 0:
        current_value += 10
    elif i % 3 == 1:
        current_value += 8
    else:
        current_value += 10
    anchor_pattern.append(current_value * 0.9)  # Apply the lift factor of 0.9

# Create a zoomed range for better visualization
zoom_range = slice(0, 100)  # Adjust the range to zoom into the first 100 points

# Plot the Pi digit difference ratios and anchor pattern
plt.figure(figsize=(12, 6))
plt.plot(pi_diff_ratios[zoom_range], label="Pi Difference Ratios", color='blue', linewidth=1.5)
plt.plot(anchor_pattern[zoom_range], label="Anchor Pattern (10, 8, 10)", linestyle='--', color='orange', linewidth=1.5)
plt.title("Comparison of Pi Digit Difference Ratios and Anchor Pattern")
plt.xlabel("Index")
plt.ylabel("Value")
plt.legend()
plt.grid()
plt.show()
