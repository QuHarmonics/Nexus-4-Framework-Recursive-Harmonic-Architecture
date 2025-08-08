import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# Simulating a recursive lattice to analyze fold growth dynamically
def generate_waveform(iteration, grid_size=8):
    x = np.linspace(0, grid_size - 1, grid_size)
    y = np.linspace(0, grid_size - 1, grid_size)
    X, Y = np.meshgrid(x, y)
    base_wave = np.sin(2 * np.pi * X / grid_size) + np.cos(2 * np.pi * Y / grid_size)
    modulator = np.sin(iteration * np.pi / 10) * np.cos(X + Y)  # Recursive modulation
    return base_wave + modulator

# Count folds (peaks and troughs) in the waveform
def count_folds(waveform):
    folds = 0
    for row in waveform:
        peaks, _ = find_peaks(row)
        troughs, _ = find_peaks(-row)  # Troughs are inverted peaks
        folds += len(peaks) + len(troughs)
    return folds

# Analyze fold growth across iterations
max_iterations = 10
grid_size = 8
fold_counts = []
waveforms = []

for i in range(1, max_iterations + 1):
    waveform = generate_waveform(i, grid_size)
    waveforms.append(waveform)
    fold_counts.append(count_folds(waveform))

# Plot fold growth
plt.figure(figsize=(10, 5))
plt.plot(range(1, max_iterations + 1), fold_counts, marker='o', label="Fold Count")
plt.axhline(y=10, color='r', linestyle='--', label="Observed Fold Limit")
plt.title("Fold Growth Across Iterations")
plt.xlabel("Iteration")
plt.ylabel("Fold Count")
plt.legend()
plt.grid()
plt.show()

# Visualize waveforms at critical iterations (1, 5, 10)
fig, axs = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={'projection': '3d'})
for idx, iteration in enumerate([1, 5, 10]):
    axs[idx].plot_surface(
        np.arange(grid_size),
        np.arange(grid_size),
        waveforms[iteration - 1],
        cmap="viridis",
        edgecolor="k"
    )
    axs[idx].set_title(f"Waveform Iteration {iteration}")
    axs[idx].set_xlabel("X Axis")
    axs[idx].set_ylabel("Y Axis")
    axs[idx].set_zlabel("Amplitude")
plt.tight_layout()
plt.show()

# Output fold counts for analysis
import pandas as pd
import ace_tools as tools

df_folds = pd.DataFrame({
    "Iteration": range(1, max_iterations + 1),
    "Fold Count": fold_counts
})
tools.display_dataframe_to_user(name="Fold Count Across Iterations", dataframe=df_folds)
