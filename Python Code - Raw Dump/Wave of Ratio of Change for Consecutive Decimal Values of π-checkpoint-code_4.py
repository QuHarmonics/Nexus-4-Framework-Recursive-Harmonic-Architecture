import numpy as np
from mpmath import mp
import plotly.graph_objects as go

# Set precision for π
mp.dps = 1000  # Set precision to 1000 decimal places

# Extract the first N decimal places of π as integers
def get_pi_decimals(num_decimals):
    pi_str = str(mp.pi)[2:]  # Get π as a string and exclude "3."
    return [int(digit) for digit in pi_str[:num_decimals]]  # Convert first num_decimals to integers

# Compute the ratio of change between consecutive digits
def compute_ratio_of_change(digits):
    changes = [abs(digits[i] - digits[i - 1]) for i in range(1, len(digits))]
    ratios = [changes[i] / changes[i - 1] if changes[i - 1] != 0 else 0 for i in range(1, len(changes))]
    return ratios

# Apply harmonic stabilization to the ratios
def harmonize_ratios(ratios, harmonic_constant):
    return [ratio / (1 + harmonic_constant * abs(ratio - harmonic_constant)) for ratio in ratios]

# Number of π decimals to process
num_decimals = 1000

# Get π decimals and compute ratios
decimals = get_pi_decimals(num_decimals)
ratios = compute_ratio_of_change(decimals)

# Apply harmonic stabilization
harmonic_constant = 0.35
harmonized_ratios = harmonize_ratios(ratios, harmonic_constant)

# Convert harmonized values to a long string
harmonized_values_string = ''.join(f'{value:.8f}' for value in harmonized_ratios)

# Print the harmonized values as a single long string
print("Harmonized Values as String:")
print(harmonized_values_string)

# Plot original and harmonized ratios
fig = go.Figure()

# Original ratios
fig.add_trace(go.Scatter(
    y=ratios[:100],  # Plot first 100 points for visualization
    mode='lines',
    name='Original Ratios'
))

# Harmonized ratios
fig.add_trace(go.Scatter(
    y=harmonized_ratios[:100],  # Plot first 100 points for visualization
    mode='lines',
    name='Harmonized Ratios'
))

# Add labels and title
fig.update_layout(
    title='Comparison of Original and Harmonized Ratios of Change',
    xaxis_title='Index',
    yaxis_title='Ratio Value',
    template='plotly_white'
)

# Show the plot
fig.show()
