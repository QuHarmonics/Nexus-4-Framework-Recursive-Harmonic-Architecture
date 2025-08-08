import numpy as np
import matplotlib.pyplot as plt

def samson_wave_to_asm(waveform, time, harmonic_constant=0.35):
    """
    Converts a waveform into ASM code at a macro level, using Samson's framework 
    for harmonic validation and completeness checks.

    Args:
        waveform (array): Input waveform amplitudes over time.
        time (array): Corresponding time values for the waveform.
        harmonic_constant (float): The harmonic constant for completeness validation.

    Returns:
        asm_code (list): The generated ASM code.
        is_harmonized (bool): Whether the ASM code fully matches the waveform.
    """
    asm_code = []  # Store the reconstructed ASM instructions
    prev_value = 0  # Track the previous amplitude
    threshold = 0.1  # Tolerance for detecting changes

    # Generate ASM code from waveform
    for i, value in enumerate(waveform):
        if i == 0:
            asm_code.append(f"PUSH {value:.2f}")
            prev_value = value
            continue

        diff = value - prev_value
        ratio = value / prev_value if prev_value != 0 else 0

        if np.abs(diff) < threshold:
            asm_code.append("NOP")
        elif diff > threshold:
            asm_code.append(f"ADD {diff:.2f}")
        elif diff < -threshold:
            asm_code.append(f"SUB {np.abs(diff):.2f}")

        prev_value = value

    asm_code.append("POP")

    # Reconstruct waveform from ASM
    reconstructed_waveform = []
    current_value = 0
    for instruction in asm_code:
        if instruction.startswith("PUSH"):
            current_value = float(instruction.split()[1])
        elif instruction.startswith("ADD"):
            current_value += float(instruction.split()[1])
        elif instruction.startswith("SUB"):
            current_value -= float(instruction.split()[1])
        reconstructed_waveform.append(current_value)

    reconstructed_waveform = np.array(reconstructed_waveform)

    # Calculate harmonic completeness
    deviation = np.abs(waveform[:len(reconstructed_waveform)] - reconstructed_waveform)
    completeness_ratio = 1 - (np.sum(deviation) / np.sum(np.abs(waveform)))
    is_harmonized = completeness_ratio >= harmonic_constant

    return asm_code, is_harmonized, reconstructed_waveform

# Generate real waveform: modulated sine wave
time = np.linspace(0, 10, 1000)  # 10 seconds at 1000 samples/second
waveform = np.sin(2 * np.pi * 1 * time) * np.sin(2 * np.pi * 0.2 * time)

# Apply the function
asm_code, is_harmonized, reconstructed_waveform = samson_wave_to_asm(waveform, time)

# Plot original and reconstructed waveforms
plt.figure(figsize=(12, 6))
plt.plot(time, waveform, label="Original Waveform", color='blue', alpha=0.7)
plt.plot(time[:len(reconstructed_waveform)], reconstructed_waveform, '--', label="Reconstructed Waveform", color='red', alpha=0.7)
plt.title("Waveform Reconstruction from ASM Code")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.show()

# Print ASM code
print("Generated ASM Code:")
for line in asm_code[:50]:  # Show first 50 lines
    print(line)

print(f"\nHarmonic Completeness: {'Harmonized' if is_harmonized else 'Incomplete'}")
