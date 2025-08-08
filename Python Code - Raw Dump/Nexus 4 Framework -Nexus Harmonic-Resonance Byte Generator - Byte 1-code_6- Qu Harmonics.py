import numpy as np
import matplotlib.pyplot as plt

# Correct first several bytes of π (after decimal)
# Source: π = 3.141592653589793238462643...
byte1 = [1, 4, 1, 5, 9, 2, 6, 5]  # Digits 1–8
byte2 = [3, 5, 8, 9, 7, 9, 3, 2]  # Digits 9–16
byte3 = [3, 8, 4, 6, 2, 6, 4, 3]  # Digits 17–24
byte4 = [3, 8, 3, 2, 7, 9, 5, 0]  # Digits 25–32
byte5 = [2, 8, 8, 4, 1, 9, 7, 1]  # Digits 33–40
byte6 = [6, 9, 3, 9, 9, 3, 7, 5]  # Digits 41–48
byte7 = [1, 0, 5, 8, 2, 0, 9, 7]  # Digits 49–56
byte8 = [4, 5, 9, 2, 3, 0, 7, 8]  # Digits 57–64

# Prepare plotting
bits = np.arange(1, 9)
bytes_list = [byte1, byte2, byte3, byte4, byte5, byte6, byte7, byte8]
colors = plt.cm.rainbow(np.linspace(0, 1, len(bytes_list)))

# Plot
plt.figure(figsize=(12, 6))
for idx, (byte, color) in enumerate(zip(bytes_list, colors), start=1):
    plt.step(bits, byte, where='mid', label=f'Byte {idx}', linewidth=2, color=color)

plt.xticks(bits, [f'Bit {i}' for i in bits])
plt.xlabel('Bit Position')
plt.ylabel('Digit Value (0–9)')
plt.title('Nexus π Byte Waveform Spectrum')
plt.grid(True)
plt.legend(ncol=4, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
