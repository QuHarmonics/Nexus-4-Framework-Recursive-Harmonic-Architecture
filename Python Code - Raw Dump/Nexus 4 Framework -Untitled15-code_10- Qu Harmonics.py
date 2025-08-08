import plotly.graph_objects as go
from mpmath import mp

# Generate the first 300 digits of Pi (excluding the initial 3) as an example
mp.dps = 301  # Set precision to 301 digits
pi_digits = [int(d) for d in str(mp.pi)[2:302]]

# Calculate the difference ratio for Pi digits (Blue Line)
pi_diff_ratios = [abs(pi_digits[i+1] - pi_digits[i]) / pi_digits[i] if pi_digits[i] != 0 else 0 for i in range(len(pi_digits)-1)]

# Create a new transformation of Pi digits (Green Line: Log-normalized differences)
pi_log_diff = [abs((pi_digits[i+1] - pi_digits[i])) / (pi_digits[i] + 1) if pi_digits[i] != 0 else 0 for i in range(len(pi_digits)-1)]

# Create the horizontal line and dot positions
dot_positions = list(range(10, len(pi_diff_ratios), 8))  # Start at 10, spaced by 8
horizontal_line_y = 0.1  # The y-coordinate for the horizontal line and dots
anchor_dots_y = [horizontal_line_y] * len(dot_positions)

# Create the figure
fig = go.Figure()

# Add the Pi difference ratios line (Blue Line)
fig.add_trace(go.Scatter(
    x=list(range(len(pi_diff_ratios))),
    y=pi_diff_ratios,
    mode='lines',
    name='Pi Difference Ratios (Blue)',
    line=dict(color='blue')
))

# Add the new Pi log difference line (Green Line)
fig.add_trace(go.Scatter(
    x=list(range(len(pi_log_diff))),
    y=pi_log_diff,
    mode='lines',
    name='Pi Log-Normalized Differences (Green)',
    line=dict(color='green')
))

# Add the horizontal line at y=0.1
fig.add_trace(go.Scatter(
    x=[0, len(pi_diff_ratios)],
    y=[horizontal_line_y, horizontal_line_y],
    mode='lines',
    name='Horizontal Line at 0.1',
    line=dict(color='orange', dash='dash')
))

# Add anchor dots at specified positions
fig.add_trace(go.Scatter(
    x=dot_positions,
    y=anchor_dots_y,
    mode='markers',
    name='Anchor Points',
    marker=dict(color='red', size=8)
))

# Customize the layout
fig.update_layout(
    title="Comparison of Pi Digit Transformations with Anchor Dots",
    xaxis_title="Index",
    yaxis_title="Value",
    legend=dict(x=0.01, y=0.99),
    template="plotly_white"
)

# Show the figure
fig.show()
