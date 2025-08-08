import numpy as np
from scipy.integrate import solve_ivp

# Lorenz-84 system
def lorenz_84(t, state):
    x, y, z = state
    dxdt = -y**2 - z**2 + 4 * x
    dydt = x * y - x * z - y + 4 * z
    dzdt = x * y + 4 * y * z - 4 * z
    return [dxdt, dydt, dzdt]

# Kuramoto-Sivashinsky equation
def kuramoto_sivashinsky(t, u, x):
    dx = x[1] - x[0]
    dudt = np.zeros_like(u)
    for i in range(len(u)):
        ip = (i + 1) % len(u)
        im = (i - 1) % len(u)
        ip2 = (i + 2) % len(u)
        im2 = (i - 2) % len(u)
        dudt[i] = - (u[ip] - 2 * u[i] + u[im]) / dx**2 - (u[ip] - u[im]) / (2 * dx) - (u[ip2] - 2 * u[i] + u[im2]) / dx**4
        # Clip dudt to prevent overflow
        dudt[i] = np.clip(dudt[i], -1e5, 1e5)
    return dudt

# Lorenz-84 system parameters
state0 = (0, 1, 1)
t_span = (0, 1)
t_eval = np.linspace(0, 1, 100)

# Kuramoto-Sivashinsky equation parameters
x = np.linspace(0, 2 * np.pi, 800)  # Increased spatial discretization
u0 = np.sin(x)
t_span_ks = (0, 0.01)
t_eval_ks = np.linspace(0, 0.01, 100)

# Simulate the Lorenz-84 system
sol = solve_ivp(lorenz_84, t_span, state0, t_eval=t_eval, method='BDF')

# Simulate the Kuramoto-Sivashinsky equation
sol_ks = solve_ivp(lambda t, u: kuramoto_sivashinsky(t, u, x), t_span_ks, u0, t_eval=t_eval_ks, method='DOP853', atol=1e-6, rtol=1e-6)

print("Lorenz-84 system final state:", sol.y[:, -1])
print("Kuramoto-Sivashinsky equation final state:", sol_ks.y[:, -1])