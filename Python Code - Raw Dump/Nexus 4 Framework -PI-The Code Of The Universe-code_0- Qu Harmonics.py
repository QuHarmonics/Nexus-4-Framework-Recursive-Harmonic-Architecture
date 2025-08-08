import numpy as np
import matplotlib.pyplot as plt
from itertools import product

# Example input data
waveform_steps = [1, 4, 1, 5, 9, 2, 6, 5]  # ASM stack waveform example
ratios = [3, 1, 4, 1, 5, 9, 2, 6]  # π digit ratios example

# Time resolution
time = np.linspace(0, len(waveform_steps) - 1, 2)

# Normalize ASM data for amplitude and frequency
stack_amplitude = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0, 1))
stack_frequency = np.interp(waveform_steps, (min(waveform_steps), max(waveform_steps)), (0.1, 2))

# Normalize π data for modulation
pi_modulator = np.interp(ratios, (min(ratios), max(ratios)), (1, 2.5))

# Resample inputs to match the fine time resolution
stack_amplitude_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_amplitude)
stack_frequency_resampled = np.interp(time, np.linspace(0, len(waveform_steps) - 1, len(waveform_steps)), stack_frequency)
pi_modulator_resampled = np.interp(time, np.linspace(0, len(ratios) - 1, len(ratios)), pi_modulator)

# Define resampled inputs for consistency
inputs_resampled = [stack_amplitude_resampled, stack_frequency_resampled, pi_modulator_resampled]
combinations = list(product(inputs_resampled, repeat=3))  # All combinations of three inputs

# Initialize the plot
plt.figure(figsize=(22, 22))

# Iterate through all combinations and plot the resulting waveforms
for i, (amplitude, frequency, modulator) in enumerate(combinations):
    # Generate the combined waveform
    combined_waveform = amplitude * np.sin(2 * np.pi * frequency * time / len(time)) * modulator
    
    # Plot the waveform
    plt.plot(time, combined_waveform, label=f"Combo {i+1}: A/F/M")

# Chart details
plt.title("All Combinations of ASM and π Modulation (27 Combinations)")
plt.xlabel("Time (Normalized)")
plt.ylabel("Waveform Amplitude")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.grid(True)
plt.tight_layout()
plt.show()
