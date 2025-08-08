# Re-import required modules after environment reset
import numpy as np
import matplotlib.pyplot as plt

# Define byte sequences
byte1 = [1, 4, 1, 5, 9, 2, 6, 5]  # Digits 1–8
byte2 = [3, 5, 8, 9, 7, 9, 3, 2]  # Digits 9–16
byte3 = [3, 8, 4, 6, 2, 6, 4, 3]  # Digits 17–24
byte4 = [3, 8, 3, 2, 7, 9, 5, 0]  # Digits 25–32
byte5 = [2, 8, 8, 4, 1, 9, 7, 1]  # Digits 33–40
byte6 = [6, 9, 3, 9, 9, 3, 7, 5]  # Digits 41–48
byte7 = [1, 0, 5, 8, 2, 0, 9, 7]  # Digits 49–56
byte8 = [4, 5, 9, 2, 3, 0, 7, 8]  # Digits 57–64

# Define bit positions for plotting
bits = np.arange(1, 9)
bits_byte2 = bits + 1  # Offset by 2 bits
bits_byte3 = bits + 2  # Offset by 3 bits

# Plotting
plt.figure(figsize=(12, 6))
plt.step(bits, byte1, where='mid', label='Byte 1', linewidth=2, color='red')
plt.step(bits_byte2, byte2, where='mid', label='Byte 2 (offset +1)', linewidth=2, color='green')
plt.step(bits_byte3, byte3, where='mid', label='Byte 3 (offset +2)', linewidth=2, color='blue')

# Labeling and display
plt.xticks(np.arange(1, 18), [f'Bit {i}' for i in range(1, 18)])
plt.xlabel('Bit Position')
plt.ylabel('Digit Value (0–9)')
plt.title('Nexus Byte Waveform with Offsets')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
