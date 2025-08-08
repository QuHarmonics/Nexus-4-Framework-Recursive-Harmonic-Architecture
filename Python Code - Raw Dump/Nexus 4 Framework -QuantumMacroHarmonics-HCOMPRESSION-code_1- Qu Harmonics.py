import numpy as np
import matplotlib.pyplot as plt# Recursive scaling with constants
def harmonic_scaling(initial_value, iterations=10):
    sequence = [initial_value]
    for n in range(1, iterations + 1):
        next_value = sequence[-1] * np.pi / ((1 + np.sqrt(5)) / 2) * 1.5 / 0.5
        sequence.append(next_value)
    return sequence

# Simulate scaling for constant 1
scaling_sequence = harmonic_scaling(1)
plt.plot(scaling_sequence, label="Scaling Sequence")
plt.title("Recursive Scaling of Constants")
plt.xlabel("Iteration")
plt.ylabel("Value")
plt.legend()
plt.grid()
plt.show()
