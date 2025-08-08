import numpy as np

def lorenz_system(x, y, z, sigma, rho, beta):
    dxdt = sigma * (y - x)
    dydt = np.clip(x * (rho - z), -100, 100) - y
    dzdt = np.clip(x * y, -100, 100) - beta * z
    return dxdt, dydt, dzdt

def recursive_harmonic_feedback(x, y, z, sigma, rho, beta, n):
    for i in range(n):
        dxdt, dydt, dzdt = lorenz_system(x, y, z, sigma, rho, beta)
        x = x + dxdt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - x) / (i + 1) + 0.01 * np.random.randn()
        y = y + dydt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - y) / (i + 1) + 0.01 * np.random.randn()
        z = z + dzdt * 0.01 * (-0.5) * np.cos(i / np.pi) + (0.5 - z) / (i + 1) + 0.01 * np.random.randn()
        x = np.clip(x, -10, 10)
        y = np.clip(y, -10, 10)
        z = np.clip(z, -10, 10)
    return x, y, z

sigma = 10
rho = 28
beta = 8/3
x0, y0, z0 = (0, 1, 1.05)
n = 1000

x_final, y_final, z_final = recursive_harmonic_feedback(x0, y0, z0, sigma, rho, beta, n)
print("Final values:", x_final, y_final, z_final)