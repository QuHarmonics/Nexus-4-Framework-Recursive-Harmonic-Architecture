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
plt.figure(figsize=(10, 6))
plt.plot(atomic_numbers_full, atomic_masses_full, marker="o", linestyle="-", label="Atomic Mass")
plt.xlabel("Atomic Number")
plt.ylabel("Atomic Mass")
plt.title("Atomic Number vs. Atomic Mass (All 118 Elements)")
plt.grid(True)
plt.legend()
plt.show()

# Plot 2: Atomic Number vs. Rate of Change in Atomic Mass
plt.figure(figsize=(10, 6))
plt.plot(atomic_numbers_full[:-1], mass_differences_full, marker="s", linestyle="--", label="Rate of Change in Mass")
plt.xlabel("Atomic Number")
plt.ylabel("Rate of Change in Atomic Mass")
plt.title("Rate of Change in Atomic Mass Between Elements (All 118 Elements)")
plt.grid(True)
plt.legend()
plt.show()

# Plot 3: Atomic Number vs. Binary Length (Grouping Pattern?)
plt.figure(figsize=(10, 6))
plt.bar(atomic_numbers_full, binary_lengths_full, color="purple", label="Binary Length of Atomic Number")
plt.xlabel("Atomic Number")
plt.ylabel("Binary Length")
plt.title("Binary Length of Atomic Number (Grouping Pattern) (All 118 Elements)")
plt.grid(axis="y")
plt.legend()
plt.show()
