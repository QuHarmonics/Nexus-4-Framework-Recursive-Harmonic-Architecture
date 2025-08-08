import numpy as np
import matplotlib.pyplot as plt

# Initialize seed waveform
waveform = [1, -1]
zeta = 0  # Represents insertion point
iterations = 100

# Recursive harmonic generation
for i in range(iterations):
    xor_value = waveform[-1] ^ waveform[-2]  # XOR of last two values
    waveform.insert(len(waveform) // 2, zeta)  # Insert zeta into the middle
    waveform.append(xor_value)  # Append XOR value as a state changer

    # Balance by inserting a reflection
    reflection = waveform[-1] * -1
    waveform.append(reflection)

# Visualize the final waveform
plt.figure(figsize=(12, 6))
plt.plot(waveform, label="Generated Waveform")
plt.axhline(0, color='gray', linestyle='--')
plt.title("Recursive Harmonic Waveform Generation")
plt.xlabel("Index")
plt.ylabel("Amplitude")
plt.legend()
plt.show()
