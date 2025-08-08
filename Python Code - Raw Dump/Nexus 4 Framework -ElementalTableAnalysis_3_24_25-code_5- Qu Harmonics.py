
# Import necessary libraries
import pandas as pd
import numpy as np
import plotly.express as px
from scipy.interpolate import interp1d

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

# Interpolate data to 128-element sample
x_interp = np.linspace(atomic_numbers_full.min(), atomic_numbers_full.max(), 128)
x_interp = np.clip(x_interp, atomic_numbers_full.min(), atomic_numbers_full.max())  # Clip interpolated values

f_mass = interp1d(atomic_numbers_full, atomic_masses_full, kind='cubic')
f_mass_diff = interp1d(atomic_numbers_full[:-1], mass_differences_full, kind='cubic')
f_binary_length = interp1d(atomic_numbers_full, binary_lengths_full, kind='cubic')

atomic_masses_interp = f_mass(x_interp)
x_interp_clipped = np.clip(x_interp[:-1], atomic_numbers_full[:-1].min(), atomic_numbers_full[:-1].max())  # Clip values for mass differences interpolation
mass_differences_interp = f_mass_diff(x_interp_clipped)
binary_lengths_interp = f_binary_length(x_interp)

# Plot 1: Atomic Number vs. Atomic Mass
fig1 = px.line(x=x_interp, y=atomic_masses_interp, title="Atomic Number vs. Atomic Mass (Interpolated 128-element sample)")
fig1.update_layout(xaxis_title="Atomic Number", yaxis_title="Atomic Mass")
fig1.show()

# Plot 2: Atomic Number vs. Rate of Change in Atomic Mass
df_mass_diff_interp = pd.DataFrame({"Atomic Number": x_interp[:-1], "Rate of Change in Mass": mass_differences_interp})
fig2 = px.line(df_mass_diff_interp, x="Atomic Number", y="Rate of Change in Mass", title="Rate of Change in Atomic Mass Between Elements (Interpolated 128-element sample)")
fig2.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig2.update_layout(xaxis_title="Atomic Number", yaxis_title="Rate of Change in Atomic Mass")
fig2.show()

# Find elements with rate of change between 1.3 and 1.9
df_mass_diff_interp['Within Range'] = df_mass_diff_interp['Rate of Change in Mass'].between(1.3, 1.9)
elements_within_range = df_mass_diff_interp[df_mass_diff_interp['Within Range']][['Atomic Number', 'Rate of Change in Mass']]
print("Elements with rate of change between 1.3 and 1.9:")
print(elements_within_range)

# Plot 3: Atomic Number vs. Binary Length (Grouping Pattern?)
df_binary_length_interp = pd.DataFrame({"Atomic Number": x_interp, "Binary Length": binary_lengths_interp})
fig3 = px.bar(df_binary_length_interp, x="Atomic Number", y="Binary Length", title="Binary Length of Atomic Number (Grouping Pattern) (Interpolated 128-element sample)")
fig3.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig3.update_layout(xaxis_title="Atomic Number", yaxis_title="Binary Length")
fig3.show()

# Count elements by binary length
binary_length_counts = df_binary_length_interp['Binary Length'].value_counts()
print("\nBinary length counts:")
print(binary_length_counts)
