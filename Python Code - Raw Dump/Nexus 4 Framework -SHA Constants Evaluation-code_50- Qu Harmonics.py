# Parameters
max_iterations_list = [5, 10, 20, 50, 100]
results = {}

for max_iterations in max_iterations_list:
    alignment_history = []
    current_waveform = initial_waveform.copy()

    for i in range(max_iterations):
        # Update waveform based on the harmonic system
        current_waveform = np.abs(np.sin(current_waveform + harmonic_constant))
        entropy = np.std(current_waveform)
        alignment_history.append(entropy)

        # Check for stabilization
        if entropy < 1e-3:  # Near-zero entropy threshold
            break

    results[max_iterations] = {
        "final_waveform": current_waveform,
        "alignment_history": alignment_history
    }

# Visualization
plt.figure(figsize=(12, 6))
for max_iterations, data in results.items():
    plt.plot(range(1, len(data["alignment_history"]) + 1), data["alignment_history"], label=f'{max_iterations} Iterations')

plt.axhline(y=harmonic_constant, color='red', linestyle='--', label='Harmonic Target (0.35)')
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.title("Entropy Trends Across Different Iteration Counts")
plt.legend()
plt.grid(True)
plt.show()
