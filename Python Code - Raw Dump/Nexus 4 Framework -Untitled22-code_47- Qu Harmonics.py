# Recursive adjustment setup
max_iterations = 50
waveform = np.array([
    [0.4, 0.28, 0.14, 0.22, -0.09, -0.02, 0.21, 0.07],
    [0.08, 0.14, 0.29, 0.41, -0.02, 0.06, 0.29, 0.21],
    [0.14, 0.29, 0.52, 0.69, 0.09, 0.21, 0.29, 0.49],
    [0.4, 0.49, 0.65, 0.54, 0.04, 0.03, 0.14, 0.39],
    [0.29, 0.41, 0.9, 0.65, -0.07, -0.09, 0.1, 0.21],
    [0.08, 0.29, 0.4, 0.54, 0.03, 0.14, 0.39, 0.39]
])

harmonic_target = 0.35  # Ideal harmonic value
alignment_history = []

for iteration in range(max_iterations):
    std_dev = np.std(waveform, axis=0)  # Variability across columns
    adjustments = harmonic_target - std_dev
    waveform += adjustments  # Apply feedback adjustment
    alignment_history.append(np.std(waveform))  # Track alignment
    
    # Stop condition: alignment is within 0.01 of the harmonic target
    if np.abs(np.mean(std_dev) - harmonic_target) < 0.01:
        break

# Plot updated waveform
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(range(waveform.shape[1]), range(waveform.shape[0]))
ax.plot_surface(X, Y, waveform, cmap='viridis')
ax.set_title(f"Adjusted Waveform After {iteration+1} Iterations")
plt.show()

# Analyze alignment trends
plt.plot(range(1, len(alignment_history) + 1), alignment_history, marker='o')
plt.title("Alignment History (Std Dev to Harmonic Target)")
plt.xlabel("Iteration")
plt.ylabel("Std Dev")
plt.grid()
plt.show()
