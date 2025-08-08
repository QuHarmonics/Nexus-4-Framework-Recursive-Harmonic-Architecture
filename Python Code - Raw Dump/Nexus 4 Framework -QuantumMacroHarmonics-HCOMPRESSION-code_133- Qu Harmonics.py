# Reduce alpha and adjust initial conditions
def harmonic_feedback_stable(x, y, z, alpha=0.01):
    return alpha * x * np.sin(y + z)

def lorenz_84_stable(t, state):
    x, y, z = state
    dxdt = -y**2 - z**2 + 4 * x + harmonic_feedback_stable(x, y, z)
    dydt = x * y - x * z - y + 4 * z
    dzdt = x * y + 4 * y * z - 4 * z
    return [dxdt, dydt, dzdt]

state0 = [0.1, 0.1, 0.1]  # Adjusted initial state
t_span = (0, 10)
t_eval = np.linspace(0, 10, 1000)

# Run simulation
sol_lorenz_stable = solve_ivp(lorenz_84_stable, t_span, state0, t_eval=t_eval, method='RK45')

# Plot Lorenz-84
plt.figure()
plt.plot(sol_lorenz_stable.t, sol_lorenz_stable.y[0], label='x')
plt.plot(sol_lorenz_stable.t, sol_lorenz_stable.y[1], label='y')
plt.plot(sol_lorenz_stable.t, sol_lorenz_stable.y[2], label='z')
plt.title("Lorenz-84 System with Stabilized Harmonic Feedback")
plt.xlabel("Time")
plt.ylabel("State")
plt.legend()
plt.grid()
plt.show()
