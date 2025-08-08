import numpy as np
from scipy.integrate import solve_ivp

def lorenz_84(t, state):
    x, y, z = state
    dxdt = -y**2 - z**2 + 4 * x
    dydt = x * y - x * z - y + 4 * z
    dzdt = x * y + 4 * y * z - 4 * z
    return dxdt, dydt, dzdt

def event(t, state):
    return np.max(np.abs(state)) - 1e10
event.terminal = True

# Lorenz-84 system parameters
state0 = (0, 0.1, 0.1)  # Different initial condition
t_span = (0, 1)  # Reduced time span

# Simulate the Lorenz-84 system
sol = solve_ivp(lorenz_84, t_span, state0, events=event, method='LSODA')

print("Lorenz-84 system final state:", sol.y[:, -1])