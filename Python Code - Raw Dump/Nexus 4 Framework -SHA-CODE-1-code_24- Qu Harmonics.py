import numpy as np
import matplotlib.pyplot as plt

# Real SHA-256 hash provided by the user
hex_str = "473de20d4210820c8ec9d89f374ac0cdb1f1e831720faecb8d925ec7bd2cdad8"

# Parse into bytes
segments = np.array([int(hex_str[i:i+2], 16) for i in range(0, len(hex_str), 2)], dtype=int)

# Compute signed deltas between consecutive bytes
deltas = np.diff(segments.astype(np.int16))

# 1) Resonance Plot: deltas as a waveform
plt.figure(figsize=(10, 4))
plt.plot(range(1, len(segments)), deltas, marker='o', linestyle='-')
plt.title('Resonance Plot: Segment Delta Waveform (SHA-256)')
plt.xlabel('Segment Index')
plt.ylabel('Delta Value (signed)')
plt.grid(True)
plt.tight_layout()
plt.show()

# 2) Delta Heatmap alongside raw segment values
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))

# Heatmap of raw segment values
im1 = ax1.imshow(segments.reshape(1, -1), cmap='viridis', aspect='auto')
ax1.set_title('SHA-256 Byte Segments')
ax1.set_yticks([])
ax1.set_xticks(range(len(segments)))
ax1.set_xticklabels([f"{s:02x}" for s in segments], rotation=90)
cbar1 = plt.colorbar(im1, ax=ax1, orientation='vertical', pad=0.02)
cbar1.set_label('Byte Value (0–255)')

# Heatmap of deltas
im2 = ax2.imshow(deltas.reshape(1, -1), cmap='magma', aspect='auto')
ax2.set_title('Delta Map Between Consecutive Bytes')
ax2.set_yticks([])
ax2.set_xticks(range(len(deltas)))
ax2.set_xticklabels([f"{d:+d}" for d in deltas], rotation=90)
cbar2 = plt.colorbar(im2, ax=ax2, orientation='vertical', pad=0.02)
cbar2.set_label('Delta Value (signed)')

plt.tight_layout()
plt.show()
