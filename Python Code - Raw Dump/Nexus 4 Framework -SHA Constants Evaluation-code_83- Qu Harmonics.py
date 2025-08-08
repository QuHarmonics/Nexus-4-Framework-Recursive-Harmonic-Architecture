import numpy as np
import matplotlib.pyplot as plt

# Parameters
A = 0.1
omega = 0.05
phi = 0
D = 0.01
tau = 100
t = np.linspace(0, 1000, 400)  # Time from 0 to 1000 units

# WSW Model for Dark Matter Influence
S = A * np.sin(omega * t + phi) + D * np.exp(-t/tau)

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(t, S, label="Dark Matter Influence on Galaxy Formation")
plt.xlabel("Time")
plt.ylabel("Influence Level")
plt.title("Modeling Dark Matter Influence Using Mark1 WSW")
plt.legend()
plt.grid(True)
plt.show()
