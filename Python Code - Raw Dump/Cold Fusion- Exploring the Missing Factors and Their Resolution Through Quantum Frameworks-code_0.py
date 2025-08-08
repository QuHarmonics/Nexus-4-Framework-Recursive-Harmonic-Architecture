import matplotlib.pyplot as plt
import numpy as np

time = np.linspace(0, 10, 100)
energy_input = np.full_like(time, 5)
energy_output = 5 + np.tanh(time - 3)
temperature = 300 + 5 * np.sin(0.5 * time)

plt.figure(figsize=(10, 6))
plt.plot(time, energy_input, label='Energy Input (Constant)', linestyle='--')
plt.plot(time, energy_output, label='Energy Output (Fusion)')
plt.plot(time, temperature, label='Temperature')
plt.xlabel('Time (s)')
plt.ylabel('Energy (J) / Temp (K)')
plt.title('Cold Fusion Simulation')
plt.legend()
plt.grid(True)
plt.show()