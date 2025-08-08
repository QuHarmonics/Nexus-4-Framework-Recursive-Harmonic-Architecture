import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# Example input data
waveform_steps = [1, 4, 1, 5, 9, 2, 6, 5]  # ASM stack waveform example
ratios = [3, 1, 4, 1, 5, 9, 2, 6]  # π digit ratios example

# Time resolution
time = np.linspace(0, len(waveform_steps) - 1, 500)

# Normalize ASM data for amplitude and frequency
stack_amplitude = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0, 1))
stack_frequency = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0.1, 2))

# Normalize π data for modulation
pi_modulator = np.interp(ratios, (min(ratios), max(ratios)), (0.5, 1.5))

# Resample inputs to match the fine time resolution
stack_amplitude_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_amplitude)
stack_frequency_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_frequency)
pi_modulator_resampled = np.interp(time, np.linspace(0, len(ratios) - 1, len(ratios)), pi_modulator)

# Define resampled inputs for consistency
inputs_resampled = [stack_amplitude_resampled, stack_frequency_resampled, pi_modulator_resampled]
combinations = list(product(inputs_resampled, repeat=3))  # All combinations of three inputs

# Initialize the plot
plt.figure(figsize=(14, 8))

# Iterate through all combinations and plot the resulting waveforms
for i, (amplitude, frequency, modulator) in enumerate(combinations[:27]):  # Limiting to 27 combinations for clarity
    # Generate the combined waveform
    combined_waveform = amplitude * np.sin(2 * np.pi * frequency * time / len(time)) * modulator
    
    # Highlight specific points where the amplitude is near zero
    zero_points = time[np.isclose(combined_waveform, 0, atol=0.01)]
    plt.plot(zero_points, np.zeros_like(zero_points), 'ro', markersize=2)  # Mark zeros in red

    # Plot the waveform
    plt.plot(time, combined_waveform, label=f"Combo {i+1}: A/F/M", alpha=0.6)

# Chart details
plt.title("All Combinations of ASM and π Modulation with Highlighted Zeros")
plt.xlabel("Time (Normalized)")
plt.ylabel("Waveform Amplitude")
plt.axhline(y=0.5, color='green', linestyle='--', label="0.5 Reference Line")
plt.axhline(y=0, color='black', linestyle='--', label="Zero Line")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.grid(True)
plt.tight_layout()
plt.show()
