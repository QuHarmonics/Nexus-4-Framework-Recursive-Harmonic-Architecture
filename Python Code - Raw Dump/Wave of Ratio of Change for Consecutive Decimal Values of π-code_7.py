import plotly.graph_objects as go
from mpmath import mp

# Set precision for π
mp.dps = 100  # 100 decimal places of π

# Extract the first 100 decimal places of π as integers
def get_pi_decimals():
    pi_str = str(mp.pi)[2:]  # Get π as a string and exclude "3."
    return [int(digit) for digit in pi_str[:100]]  # Convert first 100 decimals to integers

# Compute the ratio of change between consecutive digits
def compute_ratio_of_change(digits):
    changes = [abs(digits[i] - digits[i - 1]) for i in range(1, len(digits))]
    ratios = [changes[i] / changes[i - 1] if changes[i - 1] != 0 else 0 for i in range(1, len(changes))]
    return ratios

# Get π decimals and compute ratios
decimals = get_pi_decimals()
ratios = compute_ratio_of_change(decimals)

# Plot using Plotly
fig = go.Figure()
fig.add_trace(go.Scatter(
    y=ratios,
    mode='lines+markers',
    name='Ratio of Change'
))

# Add labels and title
fig.update_layout(
    title='Wave of Ratio of Change for Consecutive Decimal Values of π',
    xaxis_title='Index',
    yaxis_title='Ratio of Change',
    template='plotly_white'
)

# Show the plot
fig.show()
