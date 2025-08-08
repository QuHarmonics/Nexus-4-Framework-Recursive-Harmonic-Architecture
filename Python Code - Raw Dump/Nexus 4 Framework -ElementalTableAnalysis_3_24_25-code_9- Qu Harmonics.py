import pandas as pd
import numpy as np
import plotly.express as px
from scipy.interpolate import interp1d  # Import interpolation library

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

# Nyquist Oversampling
# We will oversample using linear interpolation to increase resolution between existing atomic numbers.
oversample_factor = 10  # How much to increase the resolution

# Interpolate atomic numbers and atomic masses to a finer resolution
atomic_numbers_interp = np.linspace(atomic_numbers_full.min(), atomic_numbers_full.max(), len(atomic_numbers_full) * oversample_factor)
interpolator = interp1d(atomic_numbers_full, atomic_masses_full, kind='cubic')
atomic_masses_interp = interpolator(atomic_numbers_interp)

# Calculate rate of change after interpolation
mass_diffs_interp = np.diff(atomic_masses_interp)

# Binary lengths interpolation
binary_length_interp = np.array([len(bin(int(num))[2:]) for num in atomic_numbers_interp])

# Plot 1: Atomic Number vs. Atomic Mass (with oversample)
fig1 = px.line(x=atomic_numbers_interp, y=atomic_masses_interp, title="Atomic Number vs. Atomic Mass (All 118 Elements with Oversample)")
fig1.update_layout(xaxis_title="Atomic Number", yaxis_title="Atomic Mass")
fig1.show()

# Plot 2: Atomic Number vs. Rate of Change in Atomic Mass (with oversample)
df_mass_diff_interp = pd.DataFrame({"Atomic Number": atomic_numbers_interp[:-1], "Rate of Change in Mass": mass_diffs_interp})
fig2 = px.line(df_mass_diff_interp, x="Atomic Number", y="Rate of Change in Mass", title="Rate of Change in Atomic Mass Between Elements (With Oversample)")
fig2.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig2.update_layout(xaxis_title="Atomic Number", yaxis_title="Rate of Change in Atomic Mass")
fig2.show()

# Plot 3: Atomic Number vs. Binary Length (with oversample)
df_binary_length_interp = pd.DataFrame({"Atomic Number": atomic_numbers_interp, "Binary Length": binary_length_interp})
fig3 = px.bar(df_binary_length_interp, x="Atomic Number", y="Binary Length", title="Binary Length of Atomic Number (Grouping Pattern) (All 118 Elements with Oversample)")
fig3.add_hline(y=1.61803398875, line_dash="dash", line_color="red")  # Add golden ratio line
fig3.update_layout(xaxis_title="Atomic Number", yaxis_title="Binary Length")
fig3.show()

# Count elements by binary length (oversampled)
binary_length_counts = pd.Series(binary_length_interp).value_counts()
print("\nBinary length counts (Oversampled):")
print(binary_length_counts)
