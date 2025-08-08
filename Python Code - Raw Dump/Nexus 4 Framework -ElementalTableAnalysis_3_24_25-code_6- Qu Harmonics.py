import numpy as np
import matplotlib.pyplot as plt

# Let's assume we have the atomic numbers and mass data
atomic_numbers = np.array([i for i in range(1, 119)])  # Atomic numbers from 1 to 118
atomic_masses = np.random.rand(118) * 100  # Placeholder for atomic mass data

# Calculate the difference in atomic mass between consecutive elements
mass_diffs = np.diff(atomic_masses)

# Find where gaps might exist in the binary length grouping
binary_lengths = [len(bin(x)[2:]) for x in atomic_numbers]
grouped_data = {32: [], 64: []}

for i, num in enumerate(atomic_numbers):
    if binary_lengths[i] == 5:  # Assuming 5 bits means the element is part of the 32 grouping
        grouped_data[32].append(num)
    if binary_lengths[i] == 6:  # Assuming 6 bits means part of the 64 grouping
        grouped_data[64].append(num)

# Plotting the rate of change in atomic mass and binary groupings
plt.figure(figsize=(12, 6))

# Plot the rate of change in mass
plt.subplot(1, 2, 1)
plt.plot(atomic_numbers[1:], mass_diffs, color='blue')
plt.title('Rate of Change in Atomic Mass')
plt.xlabel('Atomic Number')
plt.ylabel('Change in Atomic Mass')

# Plot binary groupings for 32 and 64
plt.subplot(1, 2, 2)
plt.bar(grouped_data[32], [0.5]*len(grouped_data[32]), color='green', label='32-bit Group')
plt.bar(grouped_data[64], [0.5]*len(grouped_data[64]), color='orange', label='64-bit Group')
plt.title('Atomic Number Grouping by Binary Length')
plt.xlabel('Atomic Number')
plt.ylabel('Binary Length Group')
plt.legend()

plt.tight_layout()
plt.show()
