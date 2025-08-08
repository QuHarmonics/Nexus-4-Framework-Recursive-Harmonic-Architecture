import numpy as np
import matplotlib.pyplot as plt

# Constants for quantum growth
INITIAL_VALUE = 1.5
GROWTH_STEP = 1.0  # Base solid number increment
FRACTIONAL_STEP = 0.5  # Quantum fractional increment
ARRAY_SIZE = 250  # Size of the H array

# Initialize the H array
def grow_wave_harmonically():
    H = np.zeros(ARRAY_SIZE)
    H[0] = INITIAL_VALUE  # Start with the initial value

    # Incrementally grow the H array
    for i in range(1, ARRAY_SIZE):
        # Decide the growth step: solid or fractional
        if i % 2 == 0:  # Example rule: Alternate solid and fractional
            H[i] = H[i - 1] + GROWTH_STEP
        else:
            H[i] = H[i - 1] + FRACTIONAL_STEP

    return H

# Visualize the resulting H array
H = grow_wave_harmonically()
plt.figure(figsize=(10, 6))
plt.plot(H, label="H(n) Harmonics", color='blue')
plt.title("Node-by-Node Quantum Wave Growth")
plt.xlabel("Index")
plt.ylabel("H(n)")
plt.legend()
plt.grid()
plt.show()

# Dump the H array values
print("H array values:", H)
