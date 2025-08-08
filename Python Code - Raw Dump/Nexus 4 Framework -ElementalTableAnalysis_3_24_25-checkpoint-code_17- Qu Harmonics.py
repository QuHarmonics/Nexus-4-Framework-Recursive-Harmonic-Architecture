import numpy as np
import pandas as pd
import plotly.graph_objects as go
from mpmath import mp

# Set precision and get the first 128 digits of Pi
mp.dps = 1300  # extra precision for 128 digits
pi_digits = str(mp.pi)[2:515]  # Extract digits after the decimal point

# Use the 1st and 4th digits as the seed (first 128 digits)
pi_digits_subseed = pi_digits[0:512]  # First 128 digits after decimal
pi_digits_array = [int(digit) for digit in pi_digits_subseed]

# Create DataFrame similar to periodic table structure
df_pi = pd.DataFrame({
    "Atomic Number": list(range(1, len(pi_digits_array) + 1)),  # Using digit position as Atomic Number
    "Element": ["Pi"] * len(pi_digits_array),  # All elements are Pi (for consistency)
    "Symbol": ["π"] * len(pi_digits_array),  # All symbols are π (Pi)
    "Pi Digit": pi_digits_array,  # The digits of Pi themselves
})

# Compute rate of change in Pi digits (difference between consecutive digits)
df_pi["Rate of Change in Pi Digit"] = df_pi["Pi Digit"].diff()

# Plot using Plotly
fig = go.Figure()

# Add scatter plot with lines connecting points
fig.add_trace(go.Scatter(
    x=df_pi["Atomic Number"],
    y=df_pi["Rate of Change in Pi Digit"],
    mode='lines+markers+text',  # Line + markers + text for symbols
    text=df_pi["Symbol"],
    textposition="top center",
    marker=dict(size=8, color='blue'),
    line=dict(color='blue', width=2),
))

# Add a horizontal line for the golden ratio
golden_ratio = 1.61803398875
fig.add_hline(y=golden_ratio, line_dash="dash", line_color="gold", annotation_text="Golden Ratio", annotation_position="bottom right")

# Customize layout
fig.update_layout(
    title="Rate of Change in Pi Digits for the First 128 Digits of Pi",
    xaxis_title="Pi Digit Position",
    yaxis_title="Rate of Change in Pi Digit",
    showlegend=False
)

# Show the plot
fig.show()
