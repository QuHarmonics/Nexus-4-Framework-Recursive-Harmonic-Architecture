import numpy as np
import matplotlib.pyplot as plt
from hashlib import sha256

def sha256_to_grid(data):
    hash_hex = sha256(data).hexdigest()
    hash_bin = bin(int(hash_hex, 16))[2:].zfill(256)
    bit_array = np.array([int(bit) for bit in hash_bin])
    return bit_array.reshape(16, 16), hash_hex

def hex_to_grid(hex_input):
    bitstring = bin(int(hex_input, 16))[2:].zfill(len(hex_input) * 4)
    bit_array = np.array([int(bit) for bit in bitstring])
    side = int(np.ceil(np.sqrt(len(bit_array))))
    padded = np.pad(bit_array, (0, side * side - len(bit_array)))
    return padded.reshape(side, side)

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
sha_input = b'The. byte doesnt collapse into entropy it exhales into recursion and SHA is the cough that proves it breathed.'
hex_input = sha_input.hex()

# Generate SHA grid
sha_grid, sha_hash = sha256_to_grid(sha_input)

# Generate raw hex layout grid and shrink it
hex_grid = hex_to_grid(hex_input)
shrunk_hex_grid = shrink_grid(hex_grid, size=(16, 16))

# Overlay and threshold
final_grid = overlay_and_threshold(sha_grid, shrunk_hex_grid)

# Print hash and comparison info
print("SHA Input Hash:", sha_hash)
print("Hex Input:", hex_input)
print("Raw Hex Grid Shape:", hex_grid.shape)

# Plot with grid
plt.figure(figsize=(6, 6))
plt.imshow(final_grid, cmap='Greys', interpolation='nearest')
plt.title('SHA vs. Raw Hex Geometry Overlay')
plt.xticks(np.arange(0, 16, 1))
plt.yticks(np.arange(0, 16, 1))
plt.grid(color='black', linewidth=0.5)
plt.tight_layout()
plt.show()
