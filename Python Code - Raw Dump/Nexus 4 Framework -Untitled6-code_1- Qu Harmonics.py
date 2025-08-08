import numpy as np
import matplotlib.pyplot as plt

N = 1000  # upper bound
is_prime = np.ones(N+1, dtype=bool)
is_prime[:2] = False

# Exclusion tracking
exclusion_map = np.zeros(N+1, dtype=int)

for p in range(2, int(N**0.5) + 1):
    if is_prime[p]:
        is_prime[p*p:N+1:p] = False
        exclusion_map[p*p:N+1:p] += 1  # Mark exclusion pressure

# Visualization
plt.figure(figsize=(12, 6))
plt.plot(np.arange(N+1), exclusion_map, label="Exclusion Pressure")
plt.scatter(np.where(is_prime)[0], [0]*sum(is_prime), c='red', s=8, label="Primes")
plt.xlabel("n")
plt.ylabel("Exclusion Level")
plt.legend()
plt.title("Exclusion Pressure vs. Prime Distribution")
plt.show()
