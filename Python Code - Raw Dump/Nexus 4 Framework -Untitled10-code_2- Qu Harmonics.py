import numpy as np
import matplotlib.pyplot as plt

theta = np.linspace(0, 2*np.pi, 100)
circle = np.array([np.cos(theta), np.sin(theta)]).T
psi_sink = 0.35 * 2 * np.pi
plt.plot(circle[:, 0], circle[:, 1], 'b-', label='S^1')
plt.plot(np.cos(psi_sink), np.sin(psi_sink), 'ro', label='ψ-sink (0.35)')
plt.arrow(0, 0, np.cos(psi_sink), np.sin(psi_sink), color='r', head_width=0.05)
plt.legend()
plt.show()