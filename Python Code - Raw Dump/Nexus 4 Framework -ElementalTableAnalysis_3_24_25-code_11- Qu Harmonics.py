import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
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
fig = px.scatter(
    df_filtered, x="Atomic Number", y="Rate of Change in Mass", text="Symbol",
    title="Rate of Change in Atomic Mass for Elements with 7-bit Length Atomic Numbers",
    labels={"Rate of Change in Mass": "Rate of Change in Mass", "Atomic Number": "Atomic Number"},
    size_max=15
)
fig.update_traces(textposition='top center')

# Add a horizontal line for the golden ratio
golden_ratio = 1.61803398875
fig.add_hline(y=golden_ratio, line_dash="dash", line_color="gold", annotation_text="Golden Ratio", annotation_position="bottom right")

fig.show()
