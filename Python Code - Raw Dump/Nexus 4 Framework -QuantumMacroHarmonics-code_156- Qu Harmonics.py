import numpy as np
import matplotlib.pyplot as plt

# Define the Riemann zeta function approximation
def zeta(s, max_iter=1000):
    if s <= 1:  # Guard against divergence for s <= 1
        raise ValueError("s must be greater than 1 for convergence")
    result = np.sum(1 / np.arange(1, max_iter + 1)**s)
    return result

# Feedback function for twin prime harmonics
def twin_prime_feedback(x, s, alpha=0.1):
    harmonic_term = alpha * x * np.cos(x)
    return harmonic_term + zeta(s)

# Initial conditions
x0 = 1  # Starting value for x
s0 = 2  # Riemann zeta function parameter (must be >1 for convergence)

# Harmonic feedback control parameters
alpha = 0.1  # Feedback scaling factor

# Simulation time
t_max = 10  # Total simulation time
dt_values = [0.01, 0.001, 0.0001]  # Different time-step sizes for testing

# Iterate over different dt values
for dt in dt_values:
    num_steps = int(t_max / dt)
    x = np.zeros(num_steps)
    s = np.zeros(num_steps)

    # Set initial conditions
    x[0] = x0
    s[0] = s0

    # Simulate the system
    for i in range(1, num_steps):
        try:
            feedback = twin_prime_feedback(x[i - 1], s[i - 1], alpha)
            x[i] = x[i - 1] + feedback * dt
            s[i] = s[i - 1]  # Keep 's' constant if desired, or modify dynamically
        except ValueError as e:
            print(f"Error at step {i} with dt={dt}: {e}")
            break

    # Plot the results for the current dt value
    plt.plot(np.linspace(0, t_max, num_steps), x, label=f"dt = {dt}")

# Finalize and display the plot
plt.legend()
plt.xlabel('Time')
plt.ylabel('x')
plt.title('Twin Prime Harmonic Feedback Control')
plt.grid()
plt.show()
