import numpy as np
import matplotlib.pyplot as plt
from hashlib import sha256

def sha256_to_grid(data):
    hash_hex = sha256(data).hexdigest()
    hash_bin = bin(int(hash_hex, 16))[2:].zfill(256)
    bit_array = np.array([int(bit) for bit in hash_bin])
    return bit_array.reshape(16, 16)

def hex_input_to_sha_grid(hex_input):
    data_bytes = bytes.fromhex(hex_input)
    return sha256_to_grid(data_bytes)

def shrink_grid(grid, size=(16, 16)):
    shrunk_grid = np.zeros(size)
    factor_x, factor_y = grid.shape[0] / size[0], grid.shape[1] / size[1]
    for x in range(size[0]):
        for y in range(size[1]):
            shrunk_grid[x, y] = np.mean(grid[int(x * factor_x):int((x + 1) * factor_x),
                                             int(y * factor_y):int((y + 1) * factor_y)])
    return shrunk_grid

def overlay_and_threshold(grid1, grid2):
    overlay = (grid1 + grid2) / 2
    midpoint = np.median(overlay)
    return overlay > midpoint

# Example inputs
sha_input = b'14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460'
hex_input = '53AB29A7AA13FED561D1A79F1A40DF5F6E03A70B8FC9411BDE93FF8B2EF81F60D2DA42B30A2EF6A0F806C361B527EAF99088B0A90C'
text_hex_input = '3134313539323635333538393739333233383436323634333338333237393530323838343139373136393339393337353130353832303937343934343539323330373831363430363238363230383939383632383033343832353334323131373036373938323134383038363531333238323330363634373039333834343630'
# Generate grids
sha_grid = sha256_to_grid(sha_input)
hex_sha_grid = hex_input_to_sha_grid(text_hex_input)
shrunk_hex_grid = shrink_grid(hex_sha_grid)

# Overlay and threshold
final_grid = overlay_and_threshold(sha_grid, shrunk_hex_grid)

# Plot with finer grid
plt.figure(figsize=(6, 6))
plt.imshow(final_grid, cmap='Greys', interpolation='nearest')
plt.title('SHA-Hex Resonance Overlay')
plt.imshow(final_grid, cmap='Greys', interpolation='nearest', extent=(0, 16, 16, 0))
plt.xticks(np.arange(0, 17, 1))
plt.yticks(np.arange(0, 17, 1))
plt.grid(which='both', color='black', linewidth=0.25)
plt.tick_params(axis='both', which='both', length=0)
plt.tight_layout()
plt.show()
