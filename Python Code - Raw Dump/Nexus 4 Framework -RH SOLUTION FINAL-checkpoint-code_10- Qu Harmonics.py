# Reinitialize data
pi_digits = [ 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4]
zeta_real_parts = [14.134725, 21.02204, 25.010857, 30.424876, 32.935062,
                   37.586178, 40.918719, 43.327073, 48.00515, 49.773832]

# Recalculate differences
pi_differences = np.diff(pi_digits)
zeta_differences = np.diff(zeta_real_parts)

# Recalculate ratios with explicit casting to float
pi_ratio_of_differences = np.divide(
    pi_differences[1:].astype(float), 
    pi_differences[:-1].astype(float), 
    out=np.zeros_like(pi_differences[1:], dtype=float), 
    where=pi_differences[:-1] != 0
)

zeta_ratio_of_differences = np.divide(
    zeta_differences[1:].astype(float), 
    zeta_differences[:-1].astype(float), 
    out=np.zeros_like(zeta_differences[1:], dtype=float), 
    where=zeta_differences[:-1] != 0
)

# Plotting the results
fig, axs = plt.subplots(2, 1, figsize=(10, 10))

# Plot the differences and ratios for π digits
axs[0].plot(pi_differences, label="Differences (π digits)", marker="o", linestyle="--", color="blue")
axs[0].plot(pi_ratio_of_differences, label="Ratio of Differences (π digits)", marker="x", linestyle="-", color="cyan")
axs[0].set_title("Differences and Ratios of Differences for π Digits")
axs[0].legend()
axs[0].grid(True)

# Plot the differences and ratios for Zeta zeros
axs[1].plot(zeta_differences, label="Differences (Zeta zeros)", marker="o", linestyle="--", color="red")
axs[1].plot(zeta_ratio_of_differences, label="Ratio of Differences (Zeta zeros)", marker="x", linestyle="-", color="orange")
axs[1].set_title("Differences and Ratios of Differences for Zeta Zeros")
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.show()
