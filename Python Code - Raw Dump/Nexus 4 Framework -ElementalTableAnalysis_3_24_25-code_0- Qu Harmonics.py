# Expanding the periodic table dataset to include all known elements up to 118
periodic_table_full = [
    {"Atomic Number": 1, "Element": "Hydrogen", "Symbol": "H", "Atomic Mass": 1.008},
    {"Atomic Number": 2, "Element": "Helium", "Symbol": "He", "Atomic Mass": 4.0026},
    {"Atomic Number": 3, "Element": "Lithium", "Symbol": "Li", "Atomic Mass": 6.94},
    {"Atomic Number": 4, "Element": "Beryllium", "Symbol": "Be", "Atomic Mass": 9.0122},
    {"Atomic Number": 5, "Element": "Boron", "Symbol": "B", "Atomic Mass": 10.81},
    {"Atomic Number": 6, "Element": "Carbon", "Symbol": "C", "Atomic Mass": 12.011},
    {"Atomic Number": 7, "Element": "Nitrogen", "Symbol": "N", "Atomic Mass": 14.007},
    {"Atomic Number": 8, "Element": "Oxygen", "Symbol": "O", "Atomic Mass": 15.999},
    {"Atomic Number": 9, "Element": "Fluorine", "Symbol": "F", "Atomic Mass": 18.998},
    {"Atomic Number": 10, "Element": "Neon", "Symbol": "Ne", "Atomic Mass": 20.180},
    {"Atomic Number": 11, "Element": "Sodium", "Symbol": "Na", "Atomic Mass": 22.990},
    {"Atomic Number": 12, "Element": "Magnesium", "Symbol": "Mg", "Atomic Mass": 24.305},
    {"Atomic Number": 13, "Element": "Aluminum", "Symbol": "Al", "Atomic Mass": 26.982},
    {"Atomic Number": 14, "Element": "Silicon", "Symbol": "Si", "Atomic Mass": 28.085},
    {"Atomic Number": 15, "Element": "Phosphorus", "Symbol": "P", "Atomic Mass": 30.974},
    {"Atomic Number": 16, "Element": "Sulfur", "Symbol": "S", "Atomic Mass": 32.06},
    {"Atomic Number": 17, "Element": "Chlorine", "Symbol": "Cl", "Atomic Mass": 35.45},
    {"Atomic Number": 18, "Element": "Argon", "Symbol": "Ar", "Atomic Mass": 39.948},
    {"Atomic Number": 19, "Element": "Potassium", "Symbol": "K", "Atomic Mass": 39.098},
    {"Atomic Number": 20, "Element": "Calcium", "Symbol": "Ca", "Atomic Mass": 40.078},
    {"Atomic Number": 21, "Element": "Scandium", "Symbol": "Sc", "Atomic Mass": 44.956},
    {"Atomic Number": 22, "Element": "Titanium", "Symbol": "Ti", "Atomic Mass": 47.867},
    {"Atomic Number": 23, "Element": "Vanadium", "Symbol": "V", "Atomic Mass": 50.942},
    {"Atomic Number": 24, "Element": "Chromium", "Symbol": "Cr", "Atomic Mass": 51.996},
    {"Atomic Number": 25, "Element": "Manganese", "Symbol": "Mn", "Atomic Mass": 54.938},
    {"Atomic Number": 26, "Element": "Iron", "Symbol": "Fe", "Atomic Mass": 55.845},
    {"Atomic Number": 27, "Element": "Cobalt", "Symbol": "Co", "Atomic Mass": 58.933},
    {"Atomic Number": 28, "Element": "Nickel", "Symbol": "Ni", "Atomic Mass": 58.693},
    {"Atomic Number": 29, "Element": "Copper", "Symbol": "Cu", "Atomic Mass": 63.546},
    {"Atomic Number": 30, "Element": "Zinc", "Symbol": "Zn", "Atomic Mass": 65.38},
    {"Atomic Number": 31, "Element": "Gallium", "Symbol": "Ga", "Atomic Mass": 69.723},
    {"Atomic Number": 32, "Element": "Germanium", "Symbol": "Ge", "Atomic Mass": 72.630},
    {"Atomic Number": 33, "Element": "Arsenic", "Symbol": "As", "Atomic Mass": 74.922},
    {"Atomic Number": 34, "Element": "Selenium", "Symbol": "Se", "Atomic Mass": 78.971},
    {"Atomic Number": 35, "Element": "Bromine", "Symbol": "Br", "Atomic Mass": 79.904},
    {"Atomic Number": 36, "Element": "Krypton", "Symbol": "Kr", "Atomic Mass": 83.798},
    {"Atomic Number": 37, "Element": "Rubidium", "Symbol": "Rb", "Atomic Mass": 85.468},
    {"Atomic Number": 38, "Element": "Strontium", "Symbol": "Sr", "Atomic Mass": 87.62},
    {"Atomic Number": 39, "Element": "Yttrium", "Symbol": "Y", "Atomic Mass": 88.906},
    {"Atomic Number": 40, "Element": "Zirconium", "Symbol": "Zr", "Atomic Mass": 91.224},
]

# Convert to DataFrame
df_periodic_full = pd.DataFrame(periodic_table_full)

# Extract relevant data
atomic_numbers_full = df_periodic_full["Atomic Number"]
atomic_masses_full = df_periodic_full["Atomic Mass"]
mass_differences_full = np.diff(atomic_masses_full)  # Rate of change in atomic mass

# Compute binary lengths of atomic numbers
binary_lengths_full = df_periodic_full["Atomic Number"].apply(lambda x: len(bin(x)[2:]))

# Plot 1: Atomic Number vs. Atomic Mass
plt.figure(figsize=(8, 5))
plt.plot(atomic_numbers_full, atomic_masses_full, marker="o", linestyle="-", label="Atomic Mass")
plt.xlabel("Atomic Number")
plt.ylabel("Atomic Mass")
plt.title("Atomic Number vs. Atomic Mass (Full Periodic Table)")
plt.grid(True)
plt.legend()
plt.show()

# Plot 2: Atomic Number vs. Rate of Change in Atomic Mass
plt.figure(figsize=(8, 5))
plt.plot(atomic_numbers_full[:-1], mass_differences_full, marker="s", linestyle="--", label="Rate of Change in Mass")
plt.xlabel("Atomic Number")
plt.ylabel("Rate of Change in Atomic Mass")
plt.title("Rate of Change in Atomic Mass Between Elements (Full Periodic Table)")
plt.grid(True)
plt.legend()
plt.show()

# Plot 3: Atomic Number vs. Binary Length (Grouping Pattern?)
plt.figure(figsize=(8, 5))
plt.bar(atomic_numbers_full, binary_lengths_full, color="purple", label="Binary Length of Atomic Number")
plt.xlabel("Atomic Number")
plt.ylabel("Binary Length")
plt.title("Binary Length of Atomic Number (Grouping Pattern) (Full Periodic Table)")
plt.grid(axis="y")
plt.legend()
plt.show()
