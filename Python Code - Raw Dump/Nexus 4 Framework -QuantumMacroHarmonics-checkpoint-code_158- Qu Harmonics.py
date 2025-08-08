import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Simplified harmonic feedback function
def harmonic_feedback(x, y, z, alpha=0.1):
    # Minimal recursive influence
    return alpha * x * np.sin(y + z)

# Lorenz-84 system with minimal intervention
def lorenz_84_simple(t, state):
    x, y, z = state
    dxdt = -y**2 - z**2 + 4 * x + harmonic_feedback(x, y, z)
    dydt = x * y - x * z - y + 4 * z
    dzdt = x * y + 4 * y * z - 4 * z
    return [dxdt, dydt, dzdt]

# Kuramoto-Sivashinsky equation
def kuramoto_sivashinsky_simple(t, u, dx):
    dudt = np.zeros_like(u)
    N = len(u)
    for i in range(N):
        ip = (i + 1) % N
        im = (i - 1) % N
        dudt[i] = - (u[ip] - 2 * u[i] + u[im]) / dx**2 - u[ip] * u[im]
    return dudt

# Simulation parameters
state0 = [1, 1, 1]  # Minimal initial condition
t_span = (0, 10)
t_eval = np.linspace(0, 10, 1000)

# Kuramoto-Sivashinsky parameters
x = np.linspace(0, 2 * np.pi, 256)
dx = x[1] - x[0]
u0 = np.sin(x)

# Run simulations
sol_lorenz = solve_ivp(lorenz_84_simple, t_span, state0, t_eval=t_eval, method='RK45', atol=1e-6, rtol=1e-6)
sol_ks = solve_ivp(lambda t, u: kuramoto_sivashinsky_simple(t, u, dx), t_span, u0, t_eval=t_eval, method='RK45', atol=1e-6, rtol=1e-6)

# Plot Lorenz-84 system
plt.figure()
plt.plot(sol_lorenz.t, sol_lorenz.y[0], label='x')
plt.plot(sol_lorenz.t, sol_lorenz.y[1], label='y')
plt.plot(sol_lorenz.t, sol_lorenz.y[2], label='z')
plt.title("Lorenz-84 System with Minimal Harmonic Feedback")
plt.xlabel("Time")
plt.ylabel("State")
plt.legend()
plt.grid()

# Plot Kuramoto-Sivashinsky system
plt.figure()
for i in range(0, len(sol_ks.t), len(sol_ks.t) // 5):
    plt.plot(x, sol_ks.y[:, i], label=f't={sol_ks.t[i]:.2f}')
plt.title("Kuramoto-Sivashinsky Equation with Minimal Feedback")
plt.xlabel("Space")
plt.ylabel("u")
plt.legend()
plt.grid()

plt.show()
