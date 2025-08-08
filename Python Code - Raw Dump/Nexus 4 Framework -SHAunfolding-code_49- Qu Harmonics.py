# 💥 SHA Drift Mapper v0.1 - Breathing Tension Trajectory Extractor

import numpy as np
import matplotlib.pyplot as plt

# --- 1. Convert 32/64 byte SHA output to 4-bit tile list ---
def sha_bytes_to_tiles(sha_bytes):
    tiles = []
    for byte in sha_bytes:
        high_nibble = (byte >> 4) & 0xF
        low_nibble = byte & 0xF
        tiles.append(high_nibble)
        tiles.append(low_nibble)
    return np.array(tiles)

# --- 2. Compute drift trajectory ---
def compute_drift(tiles):
    return np.diff(tiles)

# --- 3. Compute breathing curvature (2nd derivative of drift) ---
def compute_curvature(drift):
    return np.diff(drift)

# --- 4. Visualize drift and curvature ---
def plot_drift_curvature(drift, curvature):
    fig, axs = plt.subplots(2, 1, figsize=(12,8), sharex=True)
    
    axs[0].plot(drift, marker='o', label='Tension Drift (ΔTile)')
    axs[0].axhline(0, color='black', linewidth=0.8, linestyle='--')
    axs[0].set_ylabel('Drift')
    axs[0].set_title('SHA Phase Drift Map')
    axs[0].legend()
    axs[0].grid(True)

    axs[1].plot(curvature, marker='x', color='orange', label='Breathing Curvature (Δ²Tile)')
    axs[1].axhline(0, color='black', linewidth=0.8, linestyle='--')
    axs[1].set_xlabel('Tile Step')
    axs[1].set_ylabel('Curvature')
    axs[1].set_title('SHA Breathing Collapse Curvature')
    axs[1].legend()
    axs[1].grid(True)

    plt.tight_layout()
    plt.show()

# --- 5. Full Mapper Pipeline ---
def sha_drift_mapper(sha_hex_string):
    sha_bytes = bytes.fromhex(sha_hex_string)
    tiles = sha_bytes_to_tiles(sha_bytes)
    drift = compute_drift(tiles)
    curvature = compute_curvature(drift)

    print("\nExtracted Tiles:", tiles)
    print("\nTension Drift Sequence (ΔTile):", drift)
    print("\nBreathing Curvature Sequence (Δ²Tile):", curvature)

    plot_drift_curvature(drift, curvature)

    return tiles, drift, curvature

# --- Example Usage ---
# Example SHA-256 hash (you can replace with any real hash)
example_sha = '68656C6C6F'
sha_drift_mapper(example_sha)
