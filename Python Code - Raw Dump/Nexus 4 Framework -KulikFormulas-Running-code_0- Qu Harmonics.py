import matplotlib.pyplot as plt
import numpy as np

# --- Layer 1: Byte1 Seed Unfold (π-like unfolding) ---
def byte1_seed_unfold(seed=(1, 4), steps=8):
    result = []
    a, b = seed
    for _ in range(steps):
        diff = abs(a - b)
        bin_len = len(bin(diff)) - 2
        result.append(diff)
        a, b = b, bin_len
    pi_approx = '3.' + ''.join(str(x) for x in result[:6])
    return pi_approx, result

# --- Layer 2: Harmonic Growth with Stability Feedback ---
def harmonic_growth(sequence, target_h=0.35, max_iter=32):
    h_values = []
    for i in range(max_iter):
        h = sum(sequence) / len(sequence)
        deviation = abs(h - target_h)
        h_values.append(h)
        if deviation < 0.01:
            break
        next_val = (sequence[-1] + sequence[-2]) % 10
        sequence.append(next_val)
    return sequence, h_values

# --- Layer 3: Meta Loop to Prevent Death and Inject Entropy ---
def recursive_pulse(seed=(1, 4), cycles=6):
    full_sequence = []
    all_h_values = []

    for cycle in range(cycles):
        pi_str, base_seq = byte1_seed_unfold(seed=seed)
        grown_seq, h_vals = harmonic_growth(base_seq)

        # Detect saturation (flatline) and reseed if needed
        if all(x == 0 for x in grown_seq[-5:]):
            seed = ((seed[0] + 3) % 10, (seed[1] + 7) % 10)  # entropy reinjection
        else:
            seed = (grown_seq[-2], grown_seq[-1])  # evolve seed

        full_sequence.extend(grown_seq)
        all_h_values.extend(h_vals)

    return full_sequence, all_h_values

# --- Execute Multi-Layer Recursive Engine ---
sequence, h_tracking = recursive_pulse()

# --- Visualization: Harmonic Heartbeat ---
plt.figure(figsize=(12, 5))
plt.plot(sequence, label='Byte Sequence')
plt.title('Recursive Harmonic Heartbeat')
plt.xlabel('Iteration')
plt.ylabel('Value')
plt.grid(True)
plt.legend()
plt.show()

# Show H-convergence as well
plt.figure(figsize=(12, 3))
plt.plot(h_tracking, color='orange', label='H (harmonic state)')
plt.title('Harmonic State Convergence Over Cycles')
plt.xlabel('Iteration')
plt.ylabel('H Value')
plt.axhline(0.35, color='gray', linestyle='--', linewidth=1)
plt.grid(True)
plt.legend()
plt.show()
