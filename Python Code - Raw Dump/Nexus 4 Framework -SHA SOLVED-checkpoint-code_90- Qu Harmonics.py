from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Define Lorenz attractor equations
def lorenz(state, t, sigma=10.0, rho=28.0, beta=8.0/3.0):
    x, y, z = state
    dxdt = sigma * (y - x)
    dydt = x * (rho - z) - y
    dzdt = x * y - beta * z
    return [dxdt, dydt, dzdt]

# Initial state and time range
initial_state = [0.1, 0.1, 0.1]
time = np.linspace(0, 50, 5000)

# Simulate
states = odeint(lorenz, initial_state, time)

# Plot results
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection='3d')
ax.plot(states[:, 0], states[:, 1], states[:, 2], lw=0.5)
ax.set_title("Lorenz Attractor with Dual-State Perspective")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()
