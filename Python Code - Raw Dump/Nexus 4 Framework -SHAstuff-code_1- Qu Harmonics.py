import hashlib
import matplotlib.pyplot as plt
import numpy as np

def sha256_digest(s):
    return hashlib.sha256(s.encode('utf-8')).hexdigest()

def reverse_ascii_hex(hex_str):
    # Convert to ASCII hex (each char as hex), reverse, then return as int
    ascii_hex = ''.join(format(ord(c), '02x') for c in hex_str)
    return int(ascii_hex[::-1], 16)

# Base input string
base_input = "Hello"
base_hash = sha256_digest(base_input)
base_value = reverse_ascii_hex(base_hash)

# Generate FM-modulated SHA drift map
deltas = []
labels = []
for i in range(64):
    mod_input = f"{base_input}{i}"
    mod_hash = sha256_digest(mod_input)
    mod_value = reverse_ascii_hex(mod_hash)
    delta = abs(mod_value - base_value)
    log2_delta = np.log2(float(delta)) if delta != 0 else 0
    deltas.append(log2_delta)
    labels.append(i)

# Plot the drift as FM waveform
plt.figure(figsize=(12, 6))
plt.plot(labels, deltas, marker='o', linestyle='-', linewidth=1.5)
plt.title("SHA-FM Drift Map (log₂(delta)) for 'Hello' + nonce")
plt.xlabel("Nonce Offset")
plt.ylabel("log₂(ΔSHA)")
plt.grid(True)
plt.tight_layout()
plt.show()
