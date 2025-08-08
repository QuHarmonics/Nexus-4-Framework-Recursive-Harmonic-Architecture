import numpy as np
import matplotlib.pyplot as plt

# Function to generate primes up to a given limit using the Sieve of Eratosthenes
def generate_primes(limit):
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False  # 0 and 1 are not primes
    for n in range(2, int(np.sqrt(limit)) + 1):
        if sieve[n]:
            sieve[n*n : limit + 1 : n] = False
    return np.where(sieve)[0]

# Generate primes and calculate the ratio of distances between consecutive primes
limit = 30

primes = generate_primes(limit)
prime_gaps = np.diff(primes)  # Calculate gaps between consecutive primes
ratios = prime_gaps / primes[:-1]  # Ratio of gaps to the previous prime

# Plot the primes and their gap ratios
plt.figure(figsize=(14, 8))
plt.plot(primes[:-1], ratios, label="Prime Gap Ratios", color="blue", lw=2)
plt.axhline(np.mean(ratios), color="red", linestyle="--", label="Mean Ratio")
plt.xlabel("Prime Numbers", fontsize=14)
plt.ylabel("Gap Ratio (Gap / Previous Prime)", fontsize=14)
plt.title("Prime Gap Ratios as a Function of Prime Numbers", fontsize=16)
plt.legend(fontsize=12)
plt.grid()
plt.show()
