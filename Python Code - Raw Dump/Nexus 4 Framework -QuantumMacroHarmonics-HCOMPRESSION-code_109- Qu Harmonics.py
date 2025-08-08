import numpy as np

def chuas_circuit(x, y, z, a, b, c, d):
    dxdt = a * np.clip(y - x - c * np.abs(x + 1) + c * np.abs(x - 1), -100, 100)
    dydt = x - y + z
    dzdt = -b * np.clip(y, -100, 100) - d * np.clip(z, -100, 100)
    return dxdt, dydt, dzdt

def recursive_harmonic_feedback(x, y, z, a, b, c, d, n):
    for i in range(n):
        dxdt, dydt, dzdt = chuas_circuit(x, y, z, a, b, c, d)
        x = x + dxdt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - x) / (i + 1)
        y = y + dydt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - y) / (i + 1)
        z = z + dzdt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - z) / (i + 1)
        x = np.clip(x, -10, 10)
        y = np.clip(y, -10, 10)
        z = np.clip(z, -10, 10)
    return x, y, z

a = 9
b = 100/7
c = 1/7
d = 0.078
x0, y0, z0 = (0, 0, 0)
n = 1000

x_final, y_final, z_final = recursive_harmonic_feedback(x0, y0, z0, a, b, c, d, n)
print("Final values:", x_final, y_final, z_final)