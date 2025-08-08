import matplotlib.pyplot as plt
import numpy as np

# Function to calculate difference ratios of Pi's digits
def calculate_difference_ratios(digits):
    return [abs(digits[i] - digits[i-1]) / digits[i-1] if digits[i-1] != 0 else 0 for i in range(1, len(digits))]

# Generate Pi's first 300 digits as integers
pi_digits = [int(d) for d in "31415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679"]
pi_digits += [int(d) for d in "82148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964"]

# Calculate difference ratios
difference_ratios = calculate_difference_ratios(pi_digits[:300])

# Create the constant reference ratio
reference_ratio = [0.1111] * len(difference_ratios)

# Plot
plt.figure(figsize=(14, 6))
plt.plot(range(1, len(difference_ratios) + 1), difference_ratios, label="Difference Ratios of Pi Digits", color="blue")
plt.plot(range(1, len(difference_ratios) + 1), reference_ratio, linestyle='--', label="Constant Reference Ratio (0.1111)", color="orange")

# Highlight features
plt.xlabel("Index of Pi Decimal Pairs")
plt.ylabel("Difference Ratio")
plt.title("Comparison of Difference Ratios of Pi Digits and Constant Ratio")
plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()
