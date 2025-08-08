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

# Interpolate data to 256-element sample
x_interp_256 = np.linspace(atomic_numbers_full.min(), atomic_numbers_full.max(), 256)
x_interp_256 = np.clip(x_interp_256, atomic_numbers_full.min(), atomic_numbers_full.max())  # Clip interpolated values

f_mass_256 = interp1d(atomic_numbers_full, atomic_masses_full, kind='cubic')
f_mass_diff_256 = interp1d(atomic_numbers_full[:-1], mass_differences_full, kind='cubic')
f_binary_length_256 = interp1d(atomic_numbers_full, binary_lengths_full, kind='cubic')

atomic_masses_interp_256 = f_mass_256(x_interp_256)
x_interp_clipped_256 = np.clip(x_interp_256[:-1], atomic_numbers_full[:-1].min(), atomic_numbers_full[:-1].max())  # Clip values for mass differences interpolation
mass_differences_interp_256 = f_mass_diff_256(x_interp_clipped_256)
binary_lengths_interp_256 = f_binary_length_256(x_interp_256)

# Plot 1: Atomic Number vs. Atomic Mass
fig1 = px.line(x=x_interp_256, y=atomic_masses_interp_256, title="Atomic Number vs. Atomic Mass (Interpolated 256-element sample)")
fig1.update_layout(xaxis_title="Atomic Number", yaxis_title="Atomic Mass")
fig1.show()

# Plot 2: Atomic Number vs. Rate of Change in Atomic Mass
df_mass_diff_interp_256 = pd.DataFrame({"Atomic Number": x_interp_clipped_256, "Rate of Change in Mass": mass_differences_interp_256})
df_mass_diff_interp_256['Cumulative Difference'] = df_mass_diff_interp_256['Rate of Change in Mass'].cumsum()

fig2 = px.bar(df_mass_diff_interp_256, x="Atomic Number", y="Rate of Change in Mass", title="Rate of Change in Atomic Mass Between Elements")
fig2.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig2.update_layout(xaxis_title="Atomic Number", yaxis_title="Rate of Change in Atomic Mass")
fig2.show()

# Find elements with rate of change between 1.0 and 1.999
df_mass_diff_interp_256['Within Range'] = df_mass_diff_interp_256['Rate of Change in Mass'].between(1.0, 1.999)
elements_within_range = df_mass_diff_interp_256[df_mass_diff_interp_256['Within Range']][['Atomic Number', 'Rate of Change in Mass']]
print("Elements with rate of change between 1.0 and 1.999:")
print(elements_within_range)

# Plot 3: Atomic Number vs. Binary Length (Grouping Pattern?)
df_binary_length_interp_256 = pd.DataFrame({"Atomic Number": x_interp_256, "Binary Length": binary_lengths_interp_256})
fig3 = px.bar(df_binary_length_interp_256, x="Atomic Number", y="Binary Length", title="Binary Length of Atomic Number (Grouping Pattern) (Interpolated 256-element sample)")
fig3.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig3.update_layout(xaxis_title="Atomic Number", yaxis_title="Binary Length")
fig3.show()

# Identify elements in the 64-count group (binary length of 7)
elements_in_64_count_group = df_binary_length_interp_256[df_binary_length_interp_256['Binary Length'] == 7]

# Time-adjust elements in the 64-count group
adjusted_atomic_numbers = elements_in_64_count_group['Atomic Number'].values
adjusted_mass_differences = mass_differences_interp_256[np.where(np.in1d(x_interp_clipped_256, adjusted_atomic_numbers))]

# Ensure arrays are of the same length
min_length = min(len(adjusted_atomic_numbers), len(adjusted_mass_differences))
adjusted_atomic_numbers = adjusted_atomic_numbers[:min_length]
adjusted_mass_differences = adjusted_mass_differences[:min_length]

# Plot adjusted data
df_mass_diff_adjusted = pd.DataFrame({"Atomic Number": adjusted_atomic_numbers, "Rate of Change in Mass": adjusted_mass_differences})
fig_adjusted = px.line(df_mass_diff_adjusted, x="Atomic Number", y="Rate of Change in Mass", title="Adjusted Rate of Change in Atomic Mass for 64-count Group")
fig_adjusted.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig_adjusted.update_layout(xaxis_title="Atomic Number", yaxis_title="Rate of Change in Atomic Mass")
fig_adjusted.show()

# Count elements by binary length
binary_length_counts = df_binary_length_interp_256['Binary Length'].value_counts()
print("\nBinary length counts:")
print(binary_length_counts)