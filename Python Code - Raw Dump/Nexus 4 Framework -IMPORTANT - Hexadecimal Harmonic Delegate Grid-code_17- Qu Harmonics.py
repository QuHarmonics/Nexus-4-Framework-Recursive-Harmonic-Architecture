import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from math import floor

# Step 1: Generate SHA constants from cube roots of the first 64 primes
def generate_sha_constants():
    primes = []
    num = 2
    while len(primes) < 64:
        if all(num % i != 0 for i in range(2, int(num ** 0.5) + 1)):
            primes.append(num)
        num += 1
    K = [floor((2**32) * (p ** (1/3) % 1)) for p in primes]
    return K

# Step 2: Compute adjacent ratios (K[i+1] / K[i])
def compute_ratios(K):
    return [K[i + 1] / K[i] for i in range(len(K) - 1)]

# Step 3: Optional normalization by 3 (ratios are unchanged but symbolic)
def normalize_ratios_by_3(K):
    return [(K[i + 1] / 3) / (K[i] / 3) for i in range(len(K) - 1)]

# Step 4: Create a torque matrix for analysis
def build_torque_matrix(K, ratios):
    data = {
        'K[i]': K[:-1],
        'K[i+1]': K[1:],
        'Theta (K[i+1]/K[i])': ratios
    }
    return pd.DataFrame(data)

# Step 5: Visualize harmonic field as a polar spiral
def plot_polar_ratios(ratios):
    angles = np.linspace(0, 2 * np.pi, len(ratios), endpoint=False)
    r = np.array(ratios)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(angles, r, marker='o', linestyle='-', color='teal')
    ax.set_title("SHA Constant Torque Ratios (Polar Harmonic Field)", va='bottom')
    plt.show()

# Main Execution
if __name__ == "__main__":
    K_values = generate_sha_constants()
    ratios = compute_ratios(K_values)
    normalized_ratios = normalize_ratios_by_3(K_values)
    torque_df = build_torque_matrix(K_values, ratios)

    print(torque_df.to_string(index=False))  # Print matrix to console
    plot_polar_ratios(ratios)               # Show harmonic field
