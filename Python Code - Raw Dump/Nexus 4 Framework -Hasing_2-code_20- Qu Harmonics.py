import numpy as np
import matplotlib.pyplot as plt

# Define the Harmonic Framework function
def harmonic_framework(n, alpha=1.5, target=0.5):
    H = np.zeros(n+1)
    H[0] = target
    
    for i in range(1, n+1):
        H[i] = H[i-1] * (-0.5) * np.cos(i/np.pi) + alpha * (target - H[i-1]) / (i+1)
    
    return H

# Parameters
n = 1000  # Number of iterations
alpha = 1.5  # Amplification factor
target = 0.5  # Target value

# Run the Harmonic Framework
H = harmonic_framework(n, alpha, target)

# Visualization
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(H, label='Sequence Values')
plt.axhline(y=target, color='r', linestyle='--', label='Target Value')
plt.title('Harmonic Framework Convergence')
plt.xlabel('Iteration')
plt.ylabel('Sequence Value')
plt.legend()

plt.subplot(1, 2, 2)
plt.hist(H, bins=50, alpha=0.5, label='Sequence Values')
plt.axvline(x=target, color='r', linestyle='--', label='Target Value')
plt.title('Harmonic Framework Distribution')
plt.xlabel('Sequence Value')
plt.ylabel('Frequency')
plt.legend()

plt.tight_layout()
plt.show()