import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# Input data
waveform_steps = [1, 4, 1, 5, 9, 2, 6, 5]  # Sequence of ASM-derived values
ratios = [3, 1, 4, 1, 5, 9, 2, 6]  # π-digit-derived ratios

# Time resolution for plotting
time = np.linspace(0, len(waveform_steps) - 1, 256)

# Normalize the inputs to create amplitude, frequency, and modulation
stack_amplitude = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0, 1))
stack_frequency = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0.1, 2))
pi_modulator = np.interp(ratios, (min(ratios), max(ratios)), (0.5, 1.5))

# Resample the inputs to match the fine time resolution
stack_amplitude_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_amplitude)
stack_frequency_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_frequency)
pi_modulator_resampled = np.interp(time, np.linspace(0, len(ratios) - 1, len(ratios)), pi_modulator)

# Prepare the resampled inputs
inputs_resampled = [stack_amplitude_resampled, stack_frequency_resampled, pi_modulator_resampled]
combinations = list(product(inputs_resampled, repeat=3))  # Generate combinations of amplitude, frequency, modulation

# Create the plot
plt.figure(figsize=(14, 8))

for i, (amplitude, frequency, modulator) in enumerate(combinations[:27]):  # Limit to the first 27 combinations for clarity
    # Generate the waveform
    combined_waveform = amplitude * np.sin(2 * np.pi * frequency * time / len(time)) * modulator
    
    # Highlight zero crossings
    zero_points = time[np.isclose(combined_waveform, 0, atol=0.01)]
    plt.plot(zero_points, np.zeros_like(zero_points), 'ro', markersize=2)  # Red dots for zeros

    # Plot the waveform
    plt.plot(time, combined_waveform, label=f"Combo {i+1}", alpha=0.6)

# Add labels and references
plt.title("All Combinations of ASM and π Modulation with Highlighted Zeros")
plt.xlabel("Time (Normalized)")
plt.ylabel("Waveform Amplitude")
plt.axhline(y=0.5, color='green', linestyle='--', label="0.5 Reference Line")  # Add reference line at 0.5
plt.axhline(y=0, color='black', linestyle='--', label="Zero Line")  # Add zero line
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize='small')
plt.grid(True)
plt.tight_layout()
plt.show()
