import numpy as np
import pandas as pd
import plotly.graph_objects as go
from mendeleev import element

# Generate data for all elements (1 to 118)
elements_data = []
for i in range(1, 119):
    try:
        el = element(i)
        elements_data.append({
            "Atomic Number": el.atomic_number,
            "Element": el.name,
            "Symbol": el.symbol,
            "Atomic Mass": el.atomic_weight
        })
    except:
        pass  # In case an element doesn't have data

# Create DataFrame
df = pd.DataFrame(elements_data)

# Compute rate of change in atomic mass
df["Rate of Change in Mass"] = df["Atomic Mass"].diff()

# Filter for elements with a bit length of 7 (within the 64-grouping assumption)
df_filtered = df[df["Atomic Number"].apply(lambda x: len(bin(x)[2:]) == 7)]

# Plot using Plotly
fig = go.Figure()

# Add scatter plot with lines connecting points
fig.add_trace(go.Scatter(
    x=df_filtered["Atomic Number"],
    y=df_filtered["Rate of Change in Mass"],
    mode='lines+markers+text',  # Line + markers + text for symbols
    text=df_filtered["Symbol"],
    textposition="top center",
    marker=dict(size=8, color='blue'),
    line=dict(color='blue', width=2),
))

# Add a horizontal line for the golden ratio
golden_ratio = 1.61803398875
fig.add_hline(y=golden_ratio, line_dash="dash", line_color="gold", annotation_text="Golden Ratio", annotation_position="bottom right")

# Customize layout
fig.update_layout(
    title="Rate of Change in Atomic Mass for Elements with 7-bit Length Atomic Numbers",
    xaxis_title="Atomic Number",
    yaxis_title="Rate of Change in Mass",
    showlegend=False
)

# Show the plot
fig.show()
