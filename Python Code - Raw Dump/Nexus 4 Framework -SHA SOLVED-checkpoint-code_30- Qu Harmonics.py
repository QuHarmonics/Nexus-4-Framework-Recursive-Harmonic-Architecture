import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def generate_base_waveforms(time):
    """
    Generate base waveforms: sine, square, triangle, sawtooth.
    """
    sine_wave = np.sin(2 * np.pi * 1 * time)
    square_wave = np.sign(np.sin(2 * np.pi * 1 * time))
    triangle_wave = 2 * np.abs(2 * (time % 1) - 1) - 1
    sawtooth_wave = 2 * (time % 1) - 1
    return {
        "sine": sine_wave,
        "square": square_wave,
        "triangle": triangle_wave,
        "sawtooth": sawtooth_wave
    }

def asm_for_basic_waveforms():
    """
    Define ASM instructions for basic waveforms: sine, square, triangle, sawtooth.
    """
    return {
        "sine": ["PUSH 0", "ADD 0.1", "ADD 0.1", "SUB 0.2", "SUB 0.1", "JMP_LOOP"],
        "square": ["PUSH 1", "PUSH -1", "JMP_LOOP"],
        "triangle": ["PUSH 0", "ADD 0.1", "ADD 0.1", "SUB 0.2", "SUB 0.1", "JMP_LOOP"],
        "sawtooth": ["PUSH 0", "ADD 0.1", "ADD 0.1", "MOD_RESET", "JMP_LOOP"]
    }

def decompose_waveform(waveform, time):
    """
    Decompose a complex waveform into frequency and modulation components.
    """
    freq = np.fft.rfftfreq(len(time), d=(time[1] - time[0]))
    fft_values = np.fft.rfft(waveform)
    magnitude = np.abs(fft_values)

    return freq, magnitude

def map_wave_to_asm(waveform, wave_type):
    """
    Map a base waveform to ASM instructions.
    """
    asm_code = []
    prev_value = waveform[0]
    threshold = 0.01

    for value in waveform:
        diff = value - prev_value
        if np.abs(diff) < threshold:
            asm_code.append("NOP")
        elif diff > 0:
            asm_code.append(f"ADD {diff:.2f}")
        elif diff < 0:
            asm_code.append(f"SUB {np.abs(diff):.2f}")
        prev_value = value

    asm_code.append("JMP_LOOP")  # Loop back for periodic waves
    return wave_type, asm_code

def reconstruct_waveform_from_asm(asm_code, length):
    """
    Reconstruct a waveform from ASM instructions.
    """
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

        # Looping support
        if instruction == "JMP_LOOP" and len(reconstructed_waveform) < length:
            reconstructed_waveform = reconstructed_waveform * (length // len(reconstructed_waveform) + 1)

    return np.array(reconstructed_waveform[:length])

def wave_to_asm_full(waveform, time):
    """
    Full wave-to-ASM converter.
    """
    # Generate base waveforms and primitives
    base_waveforms = generate_base_waveforms(time)
    asm_primitives = asm_for_basic_waveforms()

    # Compare waveform to base primitives
    asm_code = []
    for wave_type, base_wave in base_waveforms.items():
        if np.allclose(waveform, base_wave, atol=0.1):  # Match to the closest waveform
            asm_code = asm_primitives[wave_type]
            break

    if not asm_code:
        asm_code = ["PUSH 0", "NOP", "POP"]  # Default for unmatched waveforms

    return asm_code

# Generate a complex waveform for testing
time = np.linspace(0, 1, 100)
waveform = np.sin(2 * np.pi * 1 * time)  # A pure sine wave

# Convert waveform to ASM
asm_code = wave_to_asm_full(waveform, time)

# Reconstruct waveform from ASM
reconstructed_waveform = reconstruct_waveform_from_asm(asm_code, len(waveform))

# Plot original and reconstructed waveforms
plt.figure(figsize=(12, 6))
plt.plot(time, waveform, label="Original Waveform", color="blue")
plt.plot(time, reconstructed_waveform, '--', label="Reconstructed Waveform", color="red")
plt.title("Waveform Reconstruction from ASM Code")
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.show()

# Print ASM code
print("Generated ASM Code:")
for line in asm_code:
    print(line)
