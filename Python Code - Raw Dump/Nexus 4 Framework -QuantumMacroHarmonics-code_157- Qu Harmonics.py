import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

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
    return dudt

# Lorenz-84 parameters
state0 = (1, 1, 1)
t_span = (0, 20)
t_eval = np.linspace(0, 20, 5000)

# Kuramoto-Sivashinsky parameters
x = np.linspace(0, 2 * np.pi, 100)  # Smaller spatial grid
u0 = np.sin(x) + 0.1 * np.random.rand(len(x))  # Add perturbations
t_span_ks = (0, 1)
t_eval_ks = np.linspace(0, 1, 1000)

# Simulate the Lorenz-84 system
sol_lorenz = solve_ivp(lorenz_84, t_span, state0, t_eval=t_eval, method='RK45')

# Simulate the Kuramoto-Sivashinsky equation
try:
    sol_ks = solve_ivp(
        lambda t, u: kuramoto_sivashinsky(t, u, x), 
        t_span_ks, 
        u0, 
        t_eval=t_eval_ks, 
        method='RK45', 
        atol=1e-6, 
        rtol=1e-6
    )
    if len(sol_ks.t) == 0 or len(sol_ks.y) == 0:
        raise ValueError("Kuramoto-Sivashinsky simulation failed.")
except Exception as e:
    print(f"Error in Kuramoto-Sivashinsky simulation: {e}")
    sol_ks = None

# Plot the Lorenz-84 system
plt.figure(figsize=(10, 5))
plt.plot(sol_lorenz.t, sol_lorenz.y[0], label='x')
plt.plot(sol_lorenz.t, sol_lorenz.y[1], label='y')
plt.plot(sol_lorenz.t, sol_lorenz.y[2], label='z')
plt.xlabel('Time')
plt.ylabel('State')
plt.title('Lorenz-84 System')
plt.legend()
plt.show()

# Plot the Kuramoto-Sivashinsky equation (if simulation succeeded)
if sol_ks:
    plt.figure(figsize=(10, 5))
    for i in range(0, len(sol_ks.t), max(1, len(sol_ks.t) // 10)):  # Ensure valid range
        plt.plot(x, sol_ks.y[:, i], label=f't={sol_ks.t[i]:.2f}')
    plt.xlabel('Space')
    plt.ylabel('u')
    plt.title('Kuramoto-Sivashinsky Equation Evolution')
    plt.legend()
    plt.show()

    # Phase-space plot
    plt.figure(figsize=(8, 6))
    plt.plot(sol_ks.y[0, :], sol_ks.y[1, :], label='Attractor')
    plt.xlabel('u1')
    plt.ylabel('u2')
    plt.title('Kuramoto-Sivashinsky Equation Attractor')
    plt.legend()
    plt.show()

# Print final states
print("Lorenz-84 system final state:", sol_lorenz.y[:, -1])
if sol_ks:
    print("Kuramoto-Sivashinsky equation final state:", sol_ks.y[:, -1])
