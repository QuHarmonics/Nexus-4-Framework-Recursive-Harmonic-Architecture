import matplotlib.pyplot as plt
import numpy as np

# Data from residue tower analysis
moduli = [30, 210, 2310, 30030]
total_pairs = [8, 48, 480, 5760]
twin_pairs = [4, 9, 20, 42]
survival_percent = [100 * t / m for t, m in zip(twin_pairs, total_pairs)]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(moduli, survival_percent, marker='o', linestyle='-', color='green')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Modulus (Mₖ)')
plt.ylabel('Twin Residue Survival (%)')
plt.title('Twin Prime Residue Survival Across Modular Towers')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.tight_layout()
plt.show()
