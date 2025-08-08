import numpy as np
import matplotlib.pyplot as plt

# Byte sequences
byte1 = [1, 4, 1, 5, 9, 2, 6, 5]  # Digits 1–8
byte2 = [3, 5, 8, 9, 7, 9, 3, 2]  # Digits 9–16
byte3 = [3, 8, 4, 6, 2, 6, 4, 3]  # Digits 17–24
byte4 = [3, 8, 3, 2, 7, 9, 5, 0]  # Digits 25–32
byte5 = [2, 8, 8, 4, 1, 9, 7, 1]  # Digits 33–40
byte6 = [6, 9, 3, 9, 9, 3, 7, 5]  # Digits 41–48
byte7 = [1, 0, 5, 8, 2, 0, 9, 7]  # Digits 49–56
byte8 = [4, 5, 9, 2, 3, 0, 7, 8]  # Digits 57–64

bits = np.arange(1, 9)

# Plot
plt.figure(figsize=(10, 4))
plt.step(bits, byte1, where='mid', label='Byte 1', linewidth=2)
#plt.step(bits, byte2, where='mid', label='Byte 2', linewidth=2)
plt.step(bits, byte3, where='mid', label='Byte 3', linewidth=2)

plt.xticks(bits, [f'Bit {i}' for i in bits])
plt.xlabel('Bit Position')
plt.ylabel('Value (Waveform Amplitude)')
plt.title('Nexus Byte Waveform: Bytes 1, 2 & 3')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
