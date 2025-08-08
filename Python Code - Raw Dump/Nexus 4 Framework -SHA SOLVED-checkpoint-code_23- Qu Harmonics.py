import matplotlib.pyplot as plt

# Define the first few digits of π and compute their differences
pi_digits = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5]  # Example digits of π
differences = [abs(pi_digits[i] - pi_digits[i - 1]) for i in range(1, len(pi_digits))]

# Define ratios of differences (current difference divided by the previous one)
ratios = [differences[i] / differences[i - 1] if differences[i - 1] != 0 else 0 for i in range(1, len(differences))]

# Plot the ratios
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(ratios) + 1), ratios, marker="o", linestyle="-", label="Ratios of Differences in π Digits")
plt.axhline(0, color="gray", linestyle="--", linewidth=0.8, label="Baseline")
plt.title("Ratios of Differences Between π Digits")
plt.xlabel("Position in Sequence (n)")
plt.ylabel("Ratio of Differences (|Difference[n]/Difference[n-1]|)")
plt.legend()
plt.grid(True)
plt.show()
