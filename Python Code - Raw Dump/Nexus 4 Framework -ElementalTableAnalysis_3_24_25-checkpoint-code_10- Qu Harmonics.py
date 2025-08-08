import numpy as np
import pandas as pd
import plotly.express as px
from mendeleev import element

# Generate the periodic table dataset using mendeleev
elements = [element(i) for i in range(1, 119)]
df = pd.DataFrame({
    "Atomic Number": [e.atomic_number for e in elements],
    "Element": [e.symbol for e in elements],
    "Atomic Mass": [e.atomic_weight for e in elements]
})

# Convert atomic numbers to binary and get the bit length
df["Binary Length"] = df["Atomic Number"].apply(lambda x: len(bin(x)[2:]))

# Filter for elements with binary length of 7 (the 64-group)
df_64_group = df[df["Binary Length"] == 7].copy()

# Calculate the rate of change in atomic mass for this group
df_64_group["Rate of Change"] = df_64_group["Atomic Mass"].diff()

# Drop NaN values resulting from the diff() operation
df_64_group.dropna(inplace=True)

# Plot the rate of change for the 64-bit length group
fig = px.line(
    df_64_group,
    x="Atomic Number",
    y="Rate of Change",
    title="Rate of Change in Atomic Mass for Elements with Binary Length 7",
    labels={"Rate of Change": "Rate of Change in Atomic Mass", "Atomic Number": "Atomic Number"}
)
fig.show()
