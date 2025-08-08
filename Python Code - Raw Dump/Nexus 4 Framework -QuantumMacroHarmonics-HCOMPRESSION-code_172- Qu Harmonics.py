import numpy as np
import matplotlib.pyplot as plt

# Helper Function for Wave Expansion
def wave_expansion(seed, base, iterations):
    expanded_wave = [seed]
    for _ in range(iterations - 1):
        next_value = expanded_wave[-1] * (base + np.sin(np.sum(expanded_wave) / len(expanded_wave)))
        expanded_wave.append(next_value)
    return expanded_wave

# Parameters
iterations = 30
bases = range(2, 11)  # Base-2 to Base-10
seed = 1.0  # Initial seed value

# Collect Results
expansions = {}
for base in bases:
    expansions[base] = wave_expansion(seed, base, iterations)

# Expected Growth
expected_growths = [seed * base ** iterations for base in bases]

# Actual Growth
actual_growths = [np.sum(expansions[base]) for base in bases]

# Plotting
plt.figure(figsize=(12, 8))
for base in bases:
    plt.plot(range(iterations), expansions[base], label=f'Base-{base}')

plt.title("Wave Expansion Across Bases", fontsize=16)
plt.xlabel("Wave Index", fontsize=12)
plt.ylabel("Expanded Value", fontsize=12)
plt.legend()
plt.show()

# Display Growth Analysis
print("Growth Analysis:")
for base, actual_growth, expected_growth in zip(bases, actual_growths, expected_growths):
    print(f"Base-{base}: Actual Growth = {actual_growth:.2f}, Expected Growth = {expected_growth:.2f}")
