import numpy as np

def logistic_map(x, r):
    return r * x * (1 - x)

def recursive_harmonic_feedback(x, r, n):
    for i in range(n):
        x = logistic_map(x, r)
        x = x * (-0.5) * np.cos(i / np.pi) + (0.5 - x) / (i + 1)
    return x

r = 3.9  # Chaotic regime
x0 = 0.1
n = 1000

x_final = recursive_harmonic_feedback(x0, r, n)
print("Final value:", x_final)