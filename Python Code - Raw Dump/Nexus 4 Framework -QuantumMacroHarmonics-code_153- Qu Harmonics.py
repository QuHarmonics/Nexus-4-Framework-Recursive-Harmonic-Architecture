import numpy as np
import matplotlib.pyplot as plt

# Define the Riemann zeta function approximation
def zeta(s, max_iter=1000):
    result = 0
    for i in range(1, max_iter + 1):
        result += 1 / i**s
    return result

# Twin Prime Harmonic Feedback Function
def twin_prime_feedback(x, s, alpha=0.1):
    harmonic_term = alpha * x * np.cos(x)
    return harmonic_term + zeta(s)

# Initial conditions
x0 = 1
s0 = 1

# Harmonic feedback control parameters
alpha = 0.1

# Simulation time
t_max = 100

# Time step values for experimentation
dt_values = [0.01, 0.001, 0.0001]

# Prepare the plot
plt.figure(figsize=(10, 6))

for dt in dt_values:
    # Initialize arrays to store results
    steps = int(t_max / dt)
    x = np.zeros(steps)
    s = np.zeros(steps)

    # Set initial conditions
    x[0] = x0
    s[0] = s0

    # Simulate the system
    for i in range(1, steps):
        x[i] = x[i - 1] + twin_prime_feedback(x[i - 1], s[i - 1], alpha) * dt
        s[i] = s[i - 1] + dt

    # Plot the results for the current dt value
    plt.plot(np.linspace(0, t_max, steps), x, label=f"dt = {dt}")

# Finalize the plot
plt.legend()
plt.xlabel('Time')
plt.ylabel('x')
plt.title('Twin Prime Harmonic Feedback Control')
plt.grid(True)
plt.show()
