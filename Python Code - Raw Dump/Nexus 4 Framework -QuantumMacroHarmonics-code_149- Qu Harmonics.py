import numpy as np
import matplotlib.pyplot as plt

def zeta(s, max_iter=1000):
    result = 0
    for i in range(1, max_iter+1):
        result += 1 / i**s
    return result

def twin_prime_feedback(x, s, alpha=0.1):
    harmonic_term = alpha * x * np.cos(x)
    return harmonic_term + zeta(s)

# Initial conditions
x0 = 1
s0 = 1

# Harmonic feedback control parameters
alpha = 0.1

# Time step
dt = 0.01

# Simulation time
t_max = 100

# Initialize arrays to store results
x = np.zeros(int(t_max/dt))
s = np.zeros(int(t_max/dt))

# Set initial conditions
x[0] = x0
s[0] = s0

# Simulate the system
for i in range(1, int(t_max/dt)):
    x[i] = x[i-1] + twin_prime_feedback(x[i-1], s[i-1], alpha) * dt
    s[i] = s[i-1] + dt

# Plot the results
plt.plot(x)
plt.xlabel('Time')
plt.ylabel('x')
plt.title('Twin Prime Harmonic Feedback Control')
plt.show()