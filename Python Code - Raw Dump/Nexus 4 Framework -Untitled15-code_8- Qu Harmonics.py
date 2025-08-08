import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets import Slider

# Generate the first 300 digits of Pi (excluding the initial 3) as an example
from mpmath import mp
mp.dps = 301  # Set precision to 301 digits
pi_digits = [int(d) for d in str(mp.pi)[2:302]]

# Calculate the difference ratio for Pi digits
pi_diff_ratios = [abs(pi_digits[i+1] - pi_digits[i]) / pi_digits[i] if pi_digits[i] != 0 else 0 for i in range(len(pi_digits)-1)]

# Create the horizontal line and dot positions
dot_positions = list(range(10, len(pi_diff_ratios), 8))  # Start at 10, spaced by 8
horizontal_line_y = 0.1  # The y-coordinate for the horizontal line and dots

# Plot the Pi digit difference ratios and anchor dots with zoom capability
fig, ax = plt.subplots(figsize=(12, 6))
plt.subplots_adjust(bottom=0.25)  # Adjust space for slider

# Plot the Pi digit difference ratios
pi_line, = ax.plot(pi_diff_ratios, label="Pi Difference Ratios", color='blue', linewidth=1.5)

# Plot the horizontal line at y=0.1
ax.axhline(y=horizontal_line_y, color='orange', linestyle='--', linewidth=1, label="Horizontal Line at 0.1")

# Add dots at specified positions
dot_scatter = ax.scatter(dot_positions, [horizontal_line_y] * len(dot_positions), color='red', s=50, label="Anchor Points")

# Customize the plot
ax.set_title("Comparison of Pi Digit Difference Ratios with Anchor Dots")
ax.set_xlabel("Index")
ax.set_ylabel("Value")
ax.legend()
ax.grid()

# Add a slider for zooming
ax_slider = plt.axes([0.2, 0.1, 0.65, 0.03], facecolor='lightgray')
slider = Slider(ax_slider, 'Zoom', 0, len(pi_diff_ratios), valinit=50, valstep=1)

# Update function for slider
def update(val):
    start = int(slider.val)
    end = start + 50  # Show 50 points at a time
    ax.set_xlim(start, end)
    fig.canvas.draw_idle()

slider.on_changed(update)

plt.show()
