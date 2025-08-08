# Create a static matplotlib plot to visualize the vector-amplitude difference and drift ratio
plt.figure(figsize=(12, 6))

# Plot Vector-Amplitude Difference
plt.plot(top_freqs_sorted['Freq (Hz)'], top_freqs_sorted['Vector-Amplitude Difference'], label='Vector-Amplitude Difference', marker='o')

# Plot Drift Ratio
plt.plot(top_freqs_sorted['Freq (Hz)'], top_freqs_sorted['Drift Ratio'], label='Drift Ratio', marker='x', linestyle='--')

plt.title('Harmonic Vector-Amplitude Dynamics from SHA Polarity Wave')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Value')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
