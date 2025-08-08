import numpy as np

def rossler_attractor(x, y, z, a, b, c):
    dxdt = -y - z
    dydt = x + a * y
    dzdt = b + z * (x - c)
    return dxdt, dydt, dzdt

def recursive_harmonic_feedback(x, y, z, a, b, c, n):
    for i in range(n):
        dxdt, dydt, dzdt = rossler_attractor(x, y, z, a, b, c)
        # Scaling and numerical stabilization
        x = x + dxdt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - x) / (i + 1)
        y = y + dydt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - y) / (i + 1)
        z = z + dzdt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - z) / (i + 1)
        # Additional stabilization
        x = np.clip(x, -10, 10)
        y = np.clip(y, -10, 10)
        z = np.clip(z, -10, 10)
    return x, y, z

a = 0.2
b = 0.2
c = 5.7
x0, y0, z0 = (0.1, 0, 0)
n = 1000

x_final, y_final, z_final = recursive_harmonic_feedback(x0, y0, z0, a, b, c, n)
print("Final values:", x_final, y_final, z_final)