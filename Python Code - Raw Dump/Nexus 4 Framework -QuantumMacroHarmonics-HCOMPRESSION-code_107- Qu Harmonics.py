import numpy as np

def henon_map(x, y, a, b):
    x_new = 1 - a * x**2 + y
    y_new = b * x
    return x_new, y_new

def recursive_harmonic_feedback(x, y, a, b, n):
    for i in range(n):
        x_new, y_new = henon_map(x, y, a, b)
        x = x_new * (-0.5) * np.cos(i / np.pi) + (0.5 - x) / (i + 1)
        y = y_new * (-0.5) * np.cos(i / np.pi) + (0.5 - y) / (i + 1)
    return x, y

a = 1.4
b = 0.3
x0, y0 = (0, 0)
n = 1000

x_final, y_final = recursive_harmonic_feedback(x0, y0, a, b, n)
print("Final values:", x_final, y_final)