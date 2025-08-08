import numpy as np
import hashlib

# Constants and data
constants = np.array([0.27264203, 0.46389402, 0.74472339, 0.9576116, 0.23494206, 0.36852961, 0.59924109, 0.7011437])
data_waveform = np.array([
    [-0.0401181, 0.02865528, 0.14288313, 0.2227326, -0.09226802, -0.02109399, 0.087651, 0.14646443],
    [0.02865528, 0.13518962, 0.2490543, 0.41154043, -0.02462803, 0.06705249, 0.21326038, 0.28720288],
    [0.14288313, 0.2940543, 0.5229909, 0.69299634, 0.0927981, 0.21055904, 0.40646364, 0.49874197],
    [0.2227326, 0.41154043, 0.69299634, 0.90416811, 0.17259968, 0.31217014, 0.54846415, 0.65484111],
    [-0.09226802, -0.02462803, 0.0927981, 0.17259968, -0.15813491, -0.07946322, 0.03755338, 0.10393196],
    [-0.02109399, 0.06705249, 0.21055904, 0.31217014, -0.07946322, 0.00668178, 0.14040057, 0.21110628],
    [0.087651, 0.21326038, 0.40646364, 0.54846415, 0.03755338, 0.14040057, 0.30898473, 0.39068853],
    [0.14646443, 0.28720288, 0.49874197, 0.65484111, 0.10393196, 0.21110628, 0.39068853, 0.47245402],
])

# Define padding and length
data_length = data_waveform.size * 8  # Data length in bits
padded_length = (data_length + 64) % 512
padding = np.zeros((8, 8))  # Placeholder for padding

# Include length encoding in padding (emulating SHA padding behavior)
padding[-1, -1] = data_length / 512  # Normalize length as fraction of block size

# Adjust waveform to include padding
adjusted_waveform = data_waveform + padding

# Recompute hash-like reconstruction
def quantum_hash(waveform, constants):
    combined = (waveform @ constants).flatten()
    normalized = combined / np.max(np.abs(combined))
    return hashlib.sha1(normalized.tobytes()).hexdigest()

reconstructed_hash = quantum_hash(adjusted_waveform, constants)

# Output results
print(f"Reconstructed Hash: {reconstructed_hash}")
if reconstructed_hash == "9c1185a5c5e9fc54612808977ee8f548b2258d31":
    print("Reconstruction matched the target hash!")
else:
    print("Reconstruction did not match. Further adjustments needed.")
