import numpy as np
import matplotlib.pyplot as plt

# Synthetic example: 8-byte segment values (in hex) for demonstration
segments = np.array([0x12, 0x34, 0x29, 0x56, 0x78, 0x65, 0x43, 0x21])
# Compute deltas between consecutive segments
deltas = np.diff(segments)

# 1) Resonance Plot: deltas as a waveform
plt.figure(figsize=(8, 4))
plt.plot(range(1, len(segments)), deltas, marker='o', linestyle='-')
plt.title('Resonance Plot: Segment Delta Waveform')
plt.xlabel('Segment Index')
plt.ylabel('Delta Value (signed)')
plt.grid(True)
plt.tight_layout()
plt.show()

# 2) SHA Digest Map Prototype: heatmap of segments and deltas
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6))

# Heatmap of raw segment values
im1 = ax1.imshow(segments.reshape(1, -1), cmap='viridis', aspect='auto')
ax1.set_title('SHA Segment Byte Values')
ax1.set_yticks([])
ax1.set_xticks(range(len(segments)))
ax1.set_xticklabels([hex(x) for x in segments])
plt.colorbar(im1, ax=ax1, orientation='vertical')

# Heatmap of deltas
im2 = ax2.imshow(deltas.reshape(1, -1), cmap='magma', aspect='auto')
ax2.set_title('Delta Map Between Segments')
ax2.set_yticks([])
ax2.set_xticks(range(len(deltas)))
ax2.set_xticklabels([f"{d:+d}" for d in deltas])
plt.colorbar(im2, ax=ax2, orientation='vertical')

plt.tight_layout()
plt.show()
