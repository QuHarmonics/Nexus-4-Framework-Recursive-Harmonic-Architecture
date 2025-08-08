import plotly.graph_objects as go
from mpmath import mp

# Set precision for π
mp.dps = 2000  # Set precision for sufficient range

# Extract the first N decimal places of π as integers
def get_pi_decimals(N):
    pi_str = str(mp.pi)[2:]  # Get π as a string and exclude "3."
    return [int(digit) for digit in pi_str[:N]]  # Convert first N decimals to integers

# Compute arithmetic patterns based on differences
def compute_arithmetic_patterns(digits):
    patterns = []
    for i in range(1, len(digits)):
        pattern = digits[i] - digits[i - 1]
        patterns.append(pattern)
    return patterns

# Compute forward and reverse patterns
N = 1000
decimals = get_pi_decimals(N)

# Forward pattern
forward_patterns = compute_arithmetic_patterns(decimals)

# Reverse pattern (from the end of the N range going backward)
reverse_patterns = compute_arithmetic_patterns(decimals[::-1])

# Plot using Plotly
fig = go.Figure()

# Forward patterns
fig.add_trace(go.Scatter(
    y=forward_patterns,
    mode='lines',
    name='Forward Arithmetic Patterns'
))

# Reverse patterns
fig.add_trace(go.Scatter(
    y=reverse_patterns,
    mode='lines',
    name='Reverse Arithmetic Patterns'
))

# Add labels and title
fig.update_layout(
    title='Overlay of Forward and Reverse Arithmetic Patterns in π',
    xaxis_title='Index',
    yaxis_title='Arithmetic Pattern Value',
    template='plotly_white'
)

# Show the plot
fig.show()
