# Harmonic Recursive Field Engine
# Synthesizing: SHA collapse, Pi-Rays, golden ratio geometry, quantum lattice, turbulence, recursive feedback

import numpy as np
import matplotlib.pyplot as plt
import hashlib
import math
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# === Config ===
NUM_POINTS = 64
ITERATIONS = 256
HARMONIC_CONSTANT = 0.35
GOLDEN_RATIO = (1 + math.sqrt(5)) / 2
SEED = "1"

# === Utilities ===
def sha256_hex(data: str) -> str:
    return hashlib.sha256(data.encode()).hexdigest()

def trailing_zeros(hex_str: str) -> int:
    binary = bin(int(hex_str, 16))[2:].zfill(256)
    return len(binary) - len(binary.rstrip('0'))

def hamming_distance(hex1: str, hex2: str) -> int:
    bin1 = bin(int(hex1, 16))[2:].zfill(256)
    bin2 = bin(int(hex2, 16))[2:].zfill(256)
    return sum(c1 != c2 for c1, c2 in zip(bin1, bin2))

# === Phase-Memory Field ===
def init_energy_gradient(n):
    return np.linspace(0, 1, n)

def recursive_redistribute(energy, iterations=10, h=HARMONIC_CONSTANT):
    states = [energy.copy()]
    for _ in range(iterations):
        forward = np.roll(energy, -1)
        backward = np.roll(energy, 1)
        energy = h * backward + (1 - h) * forward
        states.append(energy.copy())
    return states

# === Quantum Turbulence ===
def add_turbulence(state, turbulence=0.2):
    noise = np.random.uniform(-turbulence, turbulence, size=state.shape)
    return state + noise

def evolve_with_turbulence(base_state, steps=10):
    states = [base_state.copy()]
    for _ in range(steps):
        new_state = add_turbulence(states[-1])
        states.append(new_state)
    return states

# === Pi-Ray Geometry ===
def generate_triangle(center, size):
    angles = np.radians([0, 120, 240])
    return np.array([
        [center[0] + size * np.cos(a), center[1] + size * np.sin(a)] for a in angles
    ])

def subdivide_triangle(triangle, level):
    if level == 0:
        return [triangle]
    else:
        midpoints = [(triangle[i] + triangle[(i + 1) % 3]) / 2 for i in range(3)]
        new_tris = [
            np.array([triangle[0], midpoints[0], midpoints[2]]),
            np.array([midpoints[0], triangle[1], midpoints[1]]),
            np.array([midpoints[2], midpoints[1], triangle[2]]),
            np.array([midpoints[0], midpoints[1], midpoints[2]]),
        ]
        result = []
        for t in new_tris:
            result.extend(subdivide_triangle(t, level - 1))
        return result

# === SHA Harmonic Series (Hex Recursion) ===
def sha_harmonic_series_hex_recursive(seed: str, steps: int):
    hashes = []
    current = seed.encode().hex()  # Step 0: "1" → "31"
    for i in range(steps):
        digest = hashlib.sha256(current.encode()).hexdigest()
        zeros = trailing_zeros(digest)
        hashes.append((i+1, current, digest, zeros))
        current = digest  # Recursively feed in as next input
    return hashes

# === Visualization ===
def plot_energy_states(states, title):
    plt.figure(figsize=(12, 8))
    for i, state in enumerate(states):
        plt.plot(range(len(state)), state, label=f"Step {i}")
    plt.title(title)
    plt.xlabel("Points")
    plt.ylabel("Energy")
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_geometry(tris):
    plt.figure(figsize=(10, 10))
    for tri in tris:
        poly = plt.Polygon(tri, edgecolor='black', facecolor='gold', alpha=0.3)
        plt.gca().add_patch(poly)
    plt.gca().set_aspect('equal')
    plt.title("Recursive Pi-Ray Geometry")
    plt.show()

def display_sha_series(series):
    print("\nSHA Harmonic Series (Hex Recursive):")
    for i, seed, digest, zeros in series:
        print(f"{i:03} | {seed:<64} → {digest} | T-Zeros: {zeros}")

# === Run Full System ===
initial_energy = init_energy_gradient(NUM_POINTS)
redistributed = recursive_redistribute(initial_energy, ITERATIONS)
turbulent_states = evolve_with_turbulence(initial_energy, ITERATIONS)
sha_series = sha_harmonic_series_hex_recursive(SEED, ITERATIONS)
center = np.array([0, 0])
triangle = generate_triangle(center, 1)
geometry = subdivide_triangle(triangle, 4)

# === Display Results ===
plot_energy_states(redistributed, "Recursive Energy Redistribution")
plot_energy_states(turbulent_states, "Turbulent Harmonic Divergence")
plot_geometry(geometry)
display_sha_series(sha_series)
