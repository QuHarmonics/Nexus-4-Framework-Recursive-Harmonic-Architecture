# Import necessary libraries
import pandas as pd
import numpy as np
import plotly.express as px

# Expanding the periodic table dataset to include all 118 known elements
from mendeleev import element

# Generate periodic table data for the first 118 elements
periodic_table_full = []
for atomic_num in range(1, 119):
    elem = element(atomic_num)
    periodic_table_full.append({
        "Atomic Number": atomic_num,
        "Element": elem.name,
        "Symbol": elem.symbol,
        "Atomic Mass": elem.atomic_weight if elem.atomic_weight else None  # Handle unknown masses
    })

# Convert to DataFrame
df_periodic_full = pd.DataFrame(periodic_table_full).dropna()  # Remove elements with unknown masses

# Extract relevant data
atomic_numbers_full = df_periodic_full["Atomic Number"]
atomic_masses_full = df_periodic_full["Atomic Mass"]
mass_differences_full = np.diff(atomic_masses_full)  # Rate of change in atomic mass

# Compute binary lengths of atomic numbers
binary_lengths_full = df_periodic_full["Atomic Number"].apply(lambda x: len(bin(x)[2:]))

# Plot 1: Atomic Number vs. Atomic Mass
fig1 = px.line(df_periodic_full, x="Atomic Number", y="Atomic Mass", title="Atomic Number vs. Atomic Mass (All 118 Elements)")
fig1.update_layout(xaxis_title="Atomic Number", yaxis_title="Atomic Mass")
fig1.show()

# Plot 2: Atomic Number vs. Rate of Change in Atomic Mass
df_mass_diff = pd.DataFrame({"Atomic Number": atomic_numbers_full[:-1], "Rate of Change in Mass": mass_differences_full})
fig2 = px.line(df_mass_diff, x="Atomic Number", y="Rate of Change in Mass", title="Rate of Change in Atomic Mass Between Elements (All 118 Elements)")
fig2.update_layout(xaxis_title="Atomic Number", yaxis_title="Rate of Change in Atomic Mass")
fig2.show()

# Plot 3: Atomic Number vs. Binary Length (Grouping Pattern?)
df_binary_length = pd.DataFrame({"Atomic Number": atomic_numbers_full, "Binary Length": binary_lengths_full})
fig3 = px.bar(df_binary_length, x="Atomic Number", y="Binary Length", title="Binary Length of Atomic Number (Grouping Pattern) (All 118 Elements)")
fig3.update_layout(xaxis_title="Atomic Number", yaxis_title="Binary Length")
fig3.show()