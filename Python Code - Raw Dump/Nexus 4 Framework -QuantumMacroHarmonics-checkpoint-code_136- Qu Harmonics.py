import numpy as np

def ikeda_map(x, y, a, b, c, d):
    x_new = a + b * x * np.cos(c - d * np.sqrt(x**2 + y**2))
    y_new = b * y * np.sin(c - d * np.sqrt(x**2 + y**2))
    return x_new, y_new

def recursive_harmonic_feedback(x, y, a, b, c, d, n):
    for i in range(n):
        x_new, y_new = ikeda_map(x, y, a, b, c, d)
        x = x_new * (-0.5) * np.cos(i / np.pi) + (0.5 - x) / (i + 1)
        y = y_new * (-0.5) * np.cos(i / np.pi) + (0.5 - y) / (i + 1)
    return x, y

a = 1
b = 0.9
c = 0.4
d = 6
x0, y0 = (0, 0)
n = 1000

x_final, y_final = recursive_harmonic_feedback(x0, y0, a, b, c, d, n)
print("Final values:", x_final, y_final)