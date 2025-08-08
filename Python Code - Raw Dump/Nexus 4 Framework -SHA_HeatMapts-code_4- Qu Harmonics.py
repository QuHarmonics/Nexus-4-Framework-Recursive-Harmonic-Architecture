import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

# Compute flip events (1 where bit changes, 0 otherwise)
flip_matrix = np.abs(np.diff(bit_array.astype(int), axis=0))

# Heatmap: rows = bits, columns = time steps
sns.heatmap(flip_matrix.T, cmap="magma", cbar=True)
plt.xlabel("Time Step")
plt.ylabel("Bit Position")
plt.title("Bit Flip Dynamics Over Time")
plt.show()