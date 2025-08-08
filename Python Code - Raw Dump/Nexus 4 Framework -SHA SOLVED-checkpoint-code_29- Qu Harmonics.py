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

# Generate waveforms
time = np.linspace(0, 1, 100)  # One-second duration
waveforms = generate_base_waveforms(time)

# Map each waveform to ASM
asm_primitives = {}
for wave_type, waveform in waveforms.items():
    wave_type, asm_code = map_wave_to_asm(waveform, wave_type)
    asm_primitives[wave_type] = asm_code

# Display ASM mappings for each base waveform
for wave_type, asm_code in asm_primitives.items():
    print(f"ASM for {wave_type.capitalize()} Waveform:")
    for line in asm_code[:20]:  # Display first 20 lines for brevity
        print(line)
    print("\n")

# Plot all base waveforms
plt.figure(figsize=(12, 8))
for wave_type, waveform in waveforms.items():
    plt.plot(time, waveform, label=wave_type.capitalize())
plt.title("Base Waveforms")
plt.xlabel("Time")
plt.ylabel("Amplitude")
plt.legend()
plt.grid(True)
plt.show()
