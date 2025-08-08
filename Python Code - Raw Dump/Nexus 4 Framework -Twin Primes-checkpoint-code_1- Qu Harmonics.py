import matplotlib.pyplot as plt
import numpy as np
from math import gcd
from functools import reduce
import textwrap

# Define primorial moduli
moduli = [30, 210, 2310, 30030]
N = 10000  # Range limit

def survives_all_filters(n, moduli):
    for M in moduli:
        if gcd(n, M) != 1 or gcd(n + 2, M) != 1:
            return False
    return True

# Generate surviving twin candidates
survivors = [n for n in range(2, N) if survives_all_filters(n, moduli)]
num_survivors = len(survivors)

# Plot results
plt.figure(figsize=(12, 5))
plt.plot(survivors, [1] * len(survivors), 'o', markersize=4, color='orange', alpha=0.7)
plt.title("Twin Prime Echo Field (Visible Survivors After Harmonic Filters)")
plt.xlabel("n (Start of Twin Pair)")
plt.ylabel("Echo Signal (1 = survives)")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

# Display survivors
formatted = textwrap.fill(", ".join(map(str, survivors[:100])), width=100)
print(formatted)
print(f"\nTotal Survivors: {num_survivors}")