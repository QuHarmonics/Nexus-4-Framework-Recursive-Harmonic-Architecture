import plotly.graph_objects as go
from mpmath import mp

# Set precision for π
mp.dps = 10000  # Adjust the precision as needed

# Extract the first N decimal places of π
def get_pi_decimals(num_decimals):
    pi_str = str(mp.pi)[2:]  # Exclude "3."
    return [int(digit) for digit in pi_str[:num_decimals]]

# Perform arithmetic operations on π digits
def calculate_arithmetic_patterns(digits):
    pattern = []
    for i in range(1, len(digits) - 1, 2):  # Step by 2 to alternate operations
        sum_val = digits[i] + digits[i + 1]
        diff_val = digits[i + 1] - digits[i]
        pattern.append(sum_val)
        pattern.append(diff_val)
    return pattern

# Get π decimals and compute patterns
num_decimals = 1000  # Adjust the number of decimals
decimals = get_pi_decimals(num_decimals)
arithmetic_pattern = calculate_arithmetic_patterns(decimals)

# Plot the results
fig = go.Figure()
fig.add_trace(go.Scatter(
    y=arithmetic_pattern[:100],  # Plot the first 100 points
    mode='lines+markers',
    name='Arithmetic Pattern'
))

# Add labels and title
fig.update_layout(
    title='Arithmetic Patterns in Digits of π',
    xaxis_title='Index',
    yaxis_title='Pattern Value',
    template='plotly_white'
)

# Show the plot
fig.show()
