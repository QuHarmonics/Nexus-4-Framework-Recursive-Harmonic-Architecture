import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def samson_fractal_growth_animated(hash_value, outer_iterations=16, middle_iterations=8, inner_iterations=3, harmonic_constant=0.35):
    """
    Samson Fractal Growth with Animation support.
    """
    frames = []  # Store frames for animation
    harmonic_resonance = harmonic_constant

    binary_hash = list(bin(int(hash_value, 16))[2:].zfill(512))

    def recursive_feedback(value, harmonic_constant, time, feedback_strength=1):
        tuned = value * (1 + feedback_strength * harmonic_constant * time)
        return 1 if tuned >= 1 else 0

    def golden_wave_influence(bit, neighbors):
        total_influence = sum(neighbors) / len(neighbors) if neighbors else 0
        return (bit + total_influence) % 2

    def grow_fractal(data, iteration, axis_multiplier):
        new_data = []
        for i in range(len(data)):
            neighbors = [
                int(data[(i - 1) % len(data)]),
                int(data[(i + 1) % len(data)])
            ]
            influenced_bit = golden_wave_influence(int(data[i]), neighbors)
            tuned_bit = recursive_feedback(
                influenced_bit, harmonic_resonance, iteration, axis_multiplier
            )
            new_data.append(int(tuned_bit))
        return new_data

    current_data = binary_hash

    for outer_iter in range(1, outer_iterations + 1):
        for middle_iter in range(1, middle_iterations + 1):
            tuned_data = current_data.copy()

            for inner_iter in range(inner_iterations):
                for axis, axis_multiplier in zip(["x", "y", "z"], [1, 2, 3]):
                    tuned_data = grow_fractal(tuned_data, middle_iter, axis_multiplier)

            harmonic_resonance *= (1 - 0.5 * (middle_iter / middle_iterations))
            current_data = [str(int(bit)) for bit in tuned_data]

            # 🌀 Capture frame after each middle iteration
            frames.append(current_data.copy())

        if len(set(current_data)) == 1:
            harmonic_resonance = harmonic_constant * (outer_iterations - outer_iter) / outer_iterations

    return frames

# 🛠 Visualization Function
def animate_fractal(frames, interval=300):
    fig, ax = plt.subplots(figsize=(12, 4))

    def update(frame_idx):
        ax.clear()
        data = frames[frame_idx]
        ax.plot(range(len(data)), list(map(int, data)), color='cyan')
        ax.set_ylim(-0.1, 1.1)
        ax.set_title(f"Samson Fractal Growth - Frame {frame_idx + 1}/{len(frames)}")
        ax.set_xlabel('Bit Index')
        ax.set_ylabel('State (0 or 1)')
        ax.grid(True)

    ani = FuncAnimation(fig, update, frames=len(frames), interval=interval, repeat=True)
    plt.tight_layout()
    plt.show()

# 🧪 Running Everything

test_hash = "185f8db32271fe25f561a6fc938b2e264306ec304eda518007d1764826381969"
frames = samson_fractal_growth_animated(test_hash)

animate_fractal(frames)
