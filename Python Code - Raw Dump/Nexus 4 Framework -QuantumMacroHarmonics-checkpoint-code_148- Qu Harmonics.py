import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the harmonic feedback loop
def harmonic_feedback(x, y, z, alpha=0.1):
    return alpha * x * np.cos(y) + alpha * z * np.sin(x)

# Lorenz-84 system parameters
state0 = (0, 0.1, 0.1)
t_span = (0, 10)
t_eval = np.linspace(0, 10, 1000)

# Define the Lorenz-84 system's equations with harmonic feedback
def lorenz_84_harmonic(t, state):
    x, y, z = state
    dxdt = -y**2 - z**2 + 4 * x + harmonic_feedback(x, y, z)
    dydt = x * y - x * z - y + 4 * z
    dzdt = x * y + 4 * y * z - 4 * z
    return dxdt, dydt, dzdt

# Simulate the Lorenz-84 system with harmonic feedback
sol = solve_ivp(lorenz_84_harmonic, t_span, state0, t_eval=t_eval, method='LSODA')

# Plot the system's behavior
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(sol.t, sol.y[0], label='x')
plt.plot(sol.t, sol.y[1], label='y')
plt.plot(sol.t, sol.y[2], label='z')
plt.xlabel('Time')
plt.ylabel('State')
plt.title('Lorenz-84 System with Harmonic Feedback')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(sol.y[0], sol.y[1], label='x-y plane')
plt.plot(sol.y[0], sol.y[2], label='x-z plane')
plt.plot(sol.y[1], sol.y[2], label='y-z plane')
plt.xlabel('State 1')
plt.ylabel('State 2')
plt.title('Lorenz-84 System Phase Portrait with Harmonic Feedback')
plt.legend()

plt.tight_layout()
plt.show()