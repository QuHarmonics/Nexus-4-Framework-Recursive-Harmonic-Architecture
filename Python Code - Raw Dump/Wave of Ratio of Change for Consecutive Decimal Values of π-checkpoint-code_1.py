import plotly.graph_objects as go
from mpmath import mp

# Set precision for π
mp.dps = 1000  # Adjust to 1000 decimal places of π or more

# Extract the first N decimal places of π as integers
def get_pi_decimals(n):
    pi_str = str(mp.pi)[2:]  # Get π as a string and exclude "3."
    return [int(digit) for digit in pi_str[:n]]  # Convert first N decimals to integers

# Compute the ratio of change between consecutive digits
def compute_ratio_of_change(digits):
    changes = [abs(digits[i] - digits[i - 1]) for i in range(1, len(digits))]
    ratios = [changes[i] / changes[i - 1] if changes[i - 1] != 0 else 0 for i in range(1, len(changes))]
    return ratios

# Get π decimals and compute ratios
num_decimals = 1000  # Define the number of decimals to process
decimals = get_pi_decimals(num_decimals)
ratios = compute_ratio_of_change(decimals)

# Plot using Plotly
fig = go.Figure()
fig.add_trace(go.Scatter(
    x=list(range(len(ratios))),  # Use indices for the x-axis
    y=ratios,
    mode='lines+markers',
    name='Ratio of Change'
))

# Add labels and title
fig.update_layout(
    title=f'Wave of Ratio of Change for First {num_decimals} Decimal Values of π',
    xaxis_title='Index',
    yaxis_title='Ratio of Change',
    template='plotly_white'
)

# Show the plot
fig.show()
