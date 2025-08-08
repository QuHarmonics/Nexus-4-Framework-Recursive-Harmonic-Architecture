import numpy as np
import hashlib
import matplotlib.pyplot as plt

# === CONFIG ===
WAVE_POINTS = 256
STEPS = 128
AMPLITUDE = 1.0
SEED = "1"

# === WAVE FORMS ===
def sine_wave(n):
    return np.sin(np.linspace(0, 2 * np.pi, n))

def square_wave(n):
    return np.sign(np.sin(np.linspace(0, 2 * np.pi, n)))

def sawtooth_wave(n):
    return 2 * (np.linspace(0, 1, n) % 1) - 1

# === HASHING ===
def sha256_raw(input_string: str) -> str:
    return hashlib.sha256(input_string.encode()).hexdigest()

# === HARMONIC ENCODER ===
def harmonic_encode(seed, steps):
    sine = sine_wave(WAVE_POINTS)
    square = square_wave(WAVE_POINTS)
    saw = sawtooth_wave(WAVE_POINTS)

    results = []
    current = seed

    for step in range(steps):
        idx = step % WAVE_POINTS

        # Combine waveform values into a harmonized input
        harmonic = sine[idx] + square[idx] + saw[idx]
        modifier = f"{current}:{harmonic:.12f}"
        hashed = sha256_raw(modifier)

        results.append((step + 1, modifier, hashed))
        current = hashed  # recursive fold

    return results

# === DISPLAY ===
def display_wave_series(series):
    for step, mod, digest in series:
        print(f"{step:03} | {mod:<40} → {digest}")

# === EXECUTE ===
wave_sha_series = harmonic_encode(SEED, STEPS)
display_wave_series(wave_sha_series)
