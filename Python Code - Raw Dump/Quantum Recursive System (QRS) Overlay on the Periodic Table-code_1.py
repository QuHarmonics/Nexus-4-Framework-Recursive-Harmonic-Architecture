# Recreate the periodic table data for the first 9 elements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define the first 9 elements
periodic_table = [
    {"Atomic Number": 1, "Element": "Hydrogen", "Symbol": "H", "Atomic Mass": 1.008},
    {"Atomic Number": 2, "Element": "Helium", "Symbol": "He", "Atomic Mass": 4.0026},
    {"Atomic Number": 3, "Element": "Lithium", "Symbol": "Li", "Atomic Mass": 6.94},
    {"Atomic Number": 4, "Element": "Beryllium", "Symbol": "Be", "Atomic Mass": 9.0122},
    {"Atomic Number": 5, "Element": "Boron", "Symbol": "B", "Atomic Mass": 10.81},
    {"Atomic Number": 6, "Element": "Carbon", "Symbol": "C", "Atomic Mass": 12.011},
    {"Atomic Number": 7, "Element": "Nitrogen", "Symbol": "N", "Atomic Mass": 14.007},
    {"Atomic Number": 8, "Element": "Oxygen", "Symbol": "O", "Atomic Mass": 15.999},
    {"Atomic Number": 9, "Element": "Fluorine", "Symbol": "F", "Atomic Mass": 18.998},
]

# Convert to DataFrame
df_periodic = pd.DataFrame(periodic_table)

# Extract relevant data
atomic_numbers = df_periodic["Atomic Number"]
atomic_masses = df_periodic["Atomic Mass"]
mass_differences = np.diff(atomic_masses)  # Rate of change in atomic mass

# Compute binary lengths of atomic numbers
binary_lengths = df_periodic["Atomic Number"].apply(lambda x: len(bin(x)[2:]))

# Plot 1: Atomic Number vs. Atomic Mass
plt.figure(figsize=(8, 5))
plt.plot(atomic_numbers, atomic_masses, marker="o", linestyle="-", label="Atomic Mass")
plt.xlabel("Atomic Number")
plt.ylabel("Atomic Mass")
plt.title("Atomic Number vs. Atomic Mass")
plt.grid(True)
plt.legend()
plt.show()

# Plot 2: Atomic Number vs. Rate of Change in Atomic Mass
plt.figure(figsize=(8, 5))
plt.plot(atomic_numbers[:-1], mass_differences, marker="s", linestyle="--", label="Rate of Change in Mass")
plt.xlabel("Atomic Number")
plt.ylabel("Rate of Change in Atomic Mass")
plt.title("Rate of Change in Atomic Mass Between Elements")
plt.grid(True)
plt.legend()
plt.show()

# Plot 3: Atomic Number vs. Binary Length (Grouping Pattern?)
plt.figure(figsize=(8, 5))
plt.bar(atomic_numbers, binary_lengths, color="purple", label="Binary Length of Atomic Number")
plt.xlabel("Atomic Number")
plt.ylabel("Binary Length")
plt.title("Binary Length of Atomic Number (Grouping Pattern)")
plt.grid(axis="y")
plt.legend()
plt.show()
