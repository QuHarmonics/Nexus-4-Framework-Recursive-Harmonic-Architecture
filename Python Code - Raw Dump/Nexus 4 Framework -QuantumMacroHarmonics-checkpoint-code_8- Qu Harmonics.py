import numpy as np
import matplotlib.pyplot as plt

# Parameters
alpha = 1.5  # Amplification factor
target = 0.5  # Convergence target
iterations = 100  # Number of iterations
initial_value = 0.5  # H(0)

# Initialize sequence
H = [initial_value]

# Recursive formula
for n in range(1, iterations + 1):
    previous = H[-1]
    correction = alpha * (target - previous) / (n + 1)
    value = previous * (-0.5) * np.cos(n / np.pi) + correction
    H.append(value)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(H, label="H(n)")
plt.axhline(y=target, color='r', linestyle='--', label="Target (0.5)")
plt.title("Convergence of H(n) to Target (0.5)")
plt.xlabel("Iteration (n)")
plt.ylabel("H(n)")
plt.legend()
plt.grid(True)
plt.show()

# Numerical validation
errors = [abs(target - h) for h in H]
print(f"Final Error: {errors[-1]}")
