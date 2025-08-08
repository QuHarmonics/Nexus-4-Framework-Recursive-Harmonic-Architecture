# Feedback loop function for iterative refinement
def feedback_loop(matrix, correction_factor=0.1):
    feedback = correction_factor * (matrix - np.mean(matrix))
    return compress_wave(matrix + feedback)

# Initialize for feedback analysis
feedback_iterations = 10
feedback_base_waveform = final_base_waveform.copy()
feedback_perturbed_waveform = final_perturbed_waveform.copy()

# Store feedback results
feedback_base_waveforms = [feedback_base_waveform]
feedback_perturbed_waveforms = [feedback_perturbed_waveform]

for _ in range(feedback_iterations):
    feedback_base_waveform = feedback_loop(feedback_base_waveform)
    feedback_perturbed_waveform = feedback_loop(feedback_perturbed_waveform)
    feedback_base_waveforms.append(feedback_base_waveform)
    feedback_perturbed_waveforms.append(feedback_perturbed_waveform)

# Visualize feedback waveforms over iterations
fig, axes = plt.subplots(2, 5, figsize=(18, 8), subplot_kw={'projection': '3d'})
fig.suptitle("Feedback Loop Wave Morphing (Base vs Perturbed)", fontsize=16)

for i, ax in enumerate(axes.flatten()[:5]):
    X, Y = np.meshgrid(range(feedback_base_waveforms[i].shape[0]), range(feedback_base_waveforms[i].shape[1]))
    ax.plot_surface(X, Y, feedback_base_waveforms[i], cmap="plasma", edgecolor="none", alpha=0.8)
    ax.set_title(f"Base Iteration {i+1}")
    ax.set_zlim(-1, 1)

for i, ax in enumerate(axes.flatten()[5:]):
    X, Y = np.meshgrid(range(feedback_perturbed_waveforms[i].shape[0]), range(feedback_perturbed_waveforms[i].shape[1]))
    ax.plot_surface(X, Y, feedback_perturbed_waveforms[i], cmap="viridis", edgecolor="none", alpha=0.8)
    ax.set_title(f"Perturbed Iteration {i+1}")
    ax.set_zlim(-1, 1)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# Compute entropy over feedback iterations
feedback_base_entropies = [compute_entropy(w) for w in feedback_base_waveforms]
feedback_perturbed_entropies = [compute_entropy(w) for w in feedback_perturbed_waveforms]

# Plot entropy during feedback refinement
plt.figure(figsize=(8, 6))
plt.plot(range(feedback_iterations + 1), feedback_base_entropies, marker='o', label='Base Entropy (Feedback)', color='blue')
plt.plot(range(feedback_iterations + 1), feedback_perturbed_entropies, marker='o', label='Perturbed Entropy (Feedback)', color='green')
plt.title("Entropy Over Feedback Iterations")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid()
plt.show()
