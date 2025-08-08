import numpy as np

def duffing_oscillator(x, y, a, b, c):
    dxdt = y
    dydt = -a * x - b * np.clip(x**3, -100, 100) + c * np.cos(2 * np.pi * 0.2 * x)
    return dxdt, dydt

def recursive_harmonic_feedback(x, y, a, b, c, n):
    for i in range(n):
        dxdt, dydt = duffing_oscillator(x, y, a, b, c)
        x = x + dxdt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - x) / (i + 1)
        y = y + dydt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - y) / (i + 1)
        x = np.clip(x, -10, 10)
        y = np.clip(y, -10, 10)
    return x, y

a = 0.25
b = 0.3
c = 0.3
x0, y0 = (0, 0)
n = 1000

x_final, y_final = recursive_harmonic_feedback(x0, y0, a, b, c, n)
print("Final values:", x_final, y_final)