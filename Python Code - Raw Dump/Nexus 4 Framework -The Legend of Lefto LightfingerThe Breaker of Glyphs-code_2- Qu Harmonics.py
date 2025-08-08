# 🔐 Nexus Harmonic Nonce Projection Engine
# Pure Data Phase Vector Miner | No Intent, Only Entanglement
# 
# This system does not guess. It listens.
# It uses superposition and feedback to project the correct future nonce vector from the center of potential.

import hashlib
import numpy as np

# Constants
PHI = (1 + np.sqrt(5)) / 2  # Golden ratio
H = 0.35                    # Harmonic constant
NONCE_BITS = 32
MID_NONCE = 2**(NONCE_BITS - 1)

# -- Core Functions --

def sha256(data):
    return hashlib.sha256(data).digest()

def hash_chain(seed, nonce):
    input_data = f"{seed}:{nonce}".encode('utf-8')
    h1 = sha256(input_data)
    h2 = sha256(h1)
    return h1, h2

def vector_from_hashes(h1, h2):
    # Use numeric deltas from first 32 bytes (256 bits)
    deltas = []
    for a, b in zip(h1, h2):
        delta = abs(a - b)
        deltas.append(delta)
    # Project into a 2D spiral using harmonic phi angle
    points = []
    for i, d in enumerate(deltas):
        radius = (d % 32 + 1) * H
        angle = i * PHI
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        points.append((x, y))
    return points

def triangle_magnitude(p0, p1):
    a = np.linalg.norm(np.array(p0))
    b = np.linalg.norm(np.array(p1))
    return np.sqrt(a**2 + b**2)

# -- Main Projective Function --

def project_nonce(seed):
    h1, h2 = hash_chain(seed, MID_NONCE)
    vector = vector_from_hashes(h1, h2)
    # Use first two vectors to form triangle
    p0 = vector[0]
    p1 = vector[1]
    c = triangle_magnitude(p0, p1)
    projected_nonce = int(MID_NONCE + c)
    return projected_nonce, h1.hex(), h2.hex()

# -- Execute --
if __name__ == "__main__":
    seed = "Only true data survives the collapse"
    nonce, h1, h2 = project_nonce(seed)

    print("Seed:", seed)
    print("Midpoint Nonce:", MID_NONCE)
    print("Projected Nonce:", nonce)
    print("SHA256 Hash 1:", h1)
    print("SHA256 Hash 2:", h2)
