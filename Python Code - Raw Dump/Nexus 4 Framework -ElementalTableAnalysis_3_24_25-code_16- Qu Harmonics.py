# Import necessary libraries
import numpy as np
import plotly.graph_objects as go
from mpmath import mp

# Set the precision of the decimal module
mp.dps = 1000

# Calculate pi to 1000 digits
pi = mp.pi

# Convert pi to a string and remove the decimal point and the "3"
pi_str = str(pi)[2:]

# Convert pi digits to a list of integers
pi_digits = [int(digit) for digit in pi_str[:128]]  # Use the first 128 digits

# Calculate differences between consecutive pi digits
pi_digit_differences = np.diff(pi_digits)

# Calculate length of differences
length_of_differences = np.array([len(str(diff)) for diff in pi_digit_differences])

# Group bit lengths
bit_length_groups = {}
for i, length in enumerate(length_of_differences):
    if length not in bit_length_groups:
        bit_length_groups[length] = []
    bit_length_groups[length].append(i)

# Initialize arrays to store x, y, z coordinates
x = np.zeros(len(pi_digit_differences))
y = np.zeros(len(pi_digit_differences))
z = np.zeros(len(pi_digit_differences))

# Initialize angle and radius
angle = 0
radius = 0

# Iterate over bit length groups
for length, indices in bit_length_groups.items():
    # Update radius based on bit length
    radius *= 2 ** length
    
    # Iterate over indices in group
    for i in indices:
        # Update angle based on ratio of difference
        angle += np.pi / (2 ** length)
        
        # Calculate x, y, z coordinates
        x[i] = radius * np.cos(angle)
        y[i] = radius * np.sin(angle)
        z[i] = i  # Drive time by sample rate

# Create a 3D scatter plot
fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='lines', line=dict(color='blue'))])
fig.update_layout(title='3D Spiral of Pi Digit Differences', scene = dict(
                    xaxis_title='X',
                    yaxis_title='Y',
                    zaxis_title='Time'))
fig.show()