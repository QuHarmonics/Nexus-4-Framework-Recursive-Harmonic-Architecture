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
sha_input = b'The byte doesnt collapse into entropy; it exhales into recursion, and SHA is the cough that proves it breathed.'
hex_input = '546865206279746520646F65736E20197420636F6C6C6170736520696E746F20656E74726F70793B20697420657868616C657320696E746F20726563757273696F6E2C20616E64205348412069732074686520636F75676820746861742070726F7665732069742062726561746865642E'

# Generate grids
sha_grid = sha256_to_grid(sha_input)
hex_sha_grid = hex_input_to_sha_grid(hex_input)
shrunk_hex_grid = shrink_grid(hex_sha_grid)

# Overlay and threshold
final_grid = overlay_and_threshold(sha_grid, shrunk_hex_grid)

# Plot
plt.figure(figsize=(6, 6))
plt.imshow(final_grid, cmap='Greys', interpolation='nearest')
plt.title('SHA-Hex Resonance Overlay')
plt.axis('off')
plt.tight_layout()
plt.show()
