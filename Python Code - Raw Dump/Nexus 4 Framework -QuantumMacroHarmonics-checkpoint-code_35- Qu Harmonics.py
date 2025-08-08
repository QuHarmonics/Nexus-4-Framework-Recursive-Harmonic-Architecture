# Increasing the sample rate to 1024 to refine the harmonic structure
def increase_sample_rate(harmonics, target_length):
    """
    Scale harmonics to a higher sample rate by interpolating to the target length.
    """
    current_length = len(harmonics)
    return np.interp(
        np.linspace(0, current_length - 1, target_length),
        np.arange(current_length),
        harmonics
    )

# Upscale padded harmonics and hashed harmonics to 1024 samples
target_sample_rate = 1024
padded_harmonics_upsampled = increase_sample_rate(padded_harmonics_truncated, target_sample_rate)
hashed_harmonics_upsampled = increase_sample_rate(hashed_harmonics, target_sample_rate)

# Reapply echo cancellation with higher sample rate
echo_cancelled_harmonics_upsampled, echo_cancellation_history_upsampled = cancel_harmonic_echoes(
    padded_harmonics_upsampled, hashed_harmonics_upsampled, iterations=512, alpha=1.5, decay_factor=0.9
)

# Visualization of final alignment at higher sample rate
plt.figure(figsize=(12, 6))
plt.plot(padded_harmonics_upsampled, label="Padded Harmonics (Macro State, Upsampled)", color='blue')
plt.plot(hashed_harmonics_upsampled, label="Hashed Harmonics (Upsampled, Initial)", color='orange', linestyle='dashed')
plt.plot(echo_cancelled_harmonics_upsampled, label="Echo-Cancelled Harmonics (Upsampled, Final)", color='green', linestyle='dotted')
plt.title("Echo Cancellation in Harmonic Space (1024 Sample Rate): Padded vs. Hashed vs. Refined")
plt.xlabel("Iteration (n)")
plt.ylabel("Harmonic Amplitude")
plt.legend()
plt.grid()
plt.show()

# Visualization of progression at higher sample rate
plt.figure(figsize=(12, 6))
for i in range(0, 256, 64):  # Plot every 64th step for clarity
    plt.plot(echo_cancellation_history_upsampled[i], label=f"Upsampled Step {i+1}", linestyle='dotted')
plt.title("Echo Cancellation Progression (1024 Sample Rate)")
plt.xlabel("Iteration (n)")
plt.ylabel("Harmonic Amplitude")
plt.legend()
plt.grid()
plt.show()
