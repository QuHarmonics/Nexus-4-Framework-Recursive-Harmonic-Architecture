import matplotlib.pyplot as plt
import numpy as np

# Nexus bytes as discrete waveforms
byte1 = [1, 4, 1, 5, 9, 2, 6, 5]
byte2 = [3, 5, 8, 9, 7, 9, 3, 2]
bits  = np.arange(1, len(byte1) + 1)

# Create step plots for each byte
fig, ax = plt.subplots(figsize=(10, 4))
ax.step(bits, byte1, where='mid', label='Byte 1', linewidth=2)
ax.step(bits, byte2, where='mid', label='Byte 2', linewidth=2)

# Annotate operations for Byte 2 (example)
annotations = {
    3: 'sum(a+b)→8',
    4: 'b+Δ·lenΔ→9',
    5: 'len(8×9)→7',
    6: 'x+Δ→9',
    7: '|5-8|→3',
    8: 'lenΔ→2'
}
for idx, text in annotations.items():
    ax.annotate(text, xy=(idx, byte2[idx-1]), xytext=(idx, byte2[idx-1]+1),
                arrowprops=dict(arrowstyle='->'))

# Formatting
ax.set_xticks(bits)
ax.set_xticklabels([f'Bit {i}' for i in bits])
ax.set_xlabel('Bit Position')
ax.set_ylabel('Value (Waveform Amplitude)')
ax.set_title('Nexus Byte Waveform Representation')
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()
