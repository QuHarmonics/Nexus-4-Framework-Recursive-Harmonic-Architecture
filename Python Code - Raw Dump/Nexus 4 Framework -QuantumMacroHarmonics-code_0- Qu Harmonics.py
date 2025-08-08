import numpy as np
import matplotlib.pyplot as plt

# Simulate negative base recursion
def negative_base_recursion(base, iterations=100, initial_value=3):
    sequence = [initial_value]
    for n in range(1, iterations + 1):
        next_value = sequence[-1] * (-base) / (n + 1)
        sequence.append(next_value)
    return sequence

# Example: Base -2
negative_sequence = negative_base_recursion(-2)
plt.plot(negative_sequence, label="Base -2")
plt.title("Negative Base Recursion")
plt.xlabel("Iteration")
plt.ylabel("Value")
plt.legend()
plt.grid()
plt.show()
