# Reinitializing the implementation due to state reset.
def samson_fractal_growth_fixed(hash_value, outer_iterations=16, middle_iterations=8, inner_iterations=3, harmonic_constant=0.35):
    """
    Incorporates Samson's Law to grow the hash as a fractal structure, keeping resonance fixed at C = 0.35.
    Outer loop adjusts, middle loop aligns, and inner loops grow the fractal across x, y, z axes.
    """
    # Initialize variables
    unfolded_data = []
    target_resonance = harmonic_constant  # Resonance remains fixed at 0.35

    # Convert the hash into binary for base structure
    binary_hash = list(bin(int(hash_value, 16))[2:].zfill(512))
    center = len(binary_hash) // 2  # Start from the center

    def recursive_feedback(value, target_resonance, feedback_strength=1):
        """Samson's recursive reflection formula to stabilize growth."""
        return value * (1 + feedback_strength * (target_resonance - value))

    def golden_wave_influence(bit, neighbors):
        """Calculates the influence of neighbors based on golden wave principles."""
        total_influence = sum(neighbors) / len(neighbors) if neighbors else 0
        return (bit + total_influence) % 2

    def grow_fractal(data, iteration, axis_multiplier):
        """Grows the fractal along a specific axis."""
        new_data = []
        for i in range(len(data)):
            neighbors = [
                int(data[(i - 1) % len(data)]),
                int(data[(i + 1) % len(data)]),
            ]
            influenced_bit = golden_wave_influence(int(data[i]), neighbors)
            tuned_bit = recursive_feedback(
                influenced_bit, target_resonance, axis_multiplier
            )
            new_data.append(int(tuned_bit))
        return new_data

    # Outer control loop
    for outer_iter in range(1, outer_iterations + 1):
        current_data = binary_hash[:]
        print(f"Outer Iteration {outer_iter}/{outer_iterations} - Starting Fractal Growth")

        # Middle alignment loop
        for middle_iter in range(1, middle_iterations + 1):
            aligned_data = []
            print(f"  Middle Iteration {middle_iter}/{middle_iterations} - Resonance Fixed at {target_resonance:.5f}")

            # Inner tuning loops (x, y, z axes)
            for inner_iter in range(inner_iterations):
                print(f"    Inner Iteration {inner_iter + 1}/{inner_iterations} - Growing Fractal")
                for axis, axis_multiplier in zip(["x", "y", "z"], [1, 2, 3]):
                    fractal_data = grow_fractal(current_data, middle_iter, axis_multiplier)
                    aligned_data.extend(fractal_data)
                    print(f"      Axis {axis} - Fractal Data (Sample): {fractal_data[:10]}")

            # Use aligned data for the next iteration
            current_data = [str(int(bit)) for bit in aligned_data]

        unfolded_data.extend(current_data)

    return unfolded_data


# Test the Samson fractal growth with fixed resonance on the SHA-256 hash
test_hash = "185f8db32271fe25f561a6fc938b2e264306ec304eda518007d1764826381969"
fractal_output_fixed = samson_fractal_growth_fixed(test_hash)

# Convert fractal output to a continuous binary string
fractal_binary_string_fixed = ''.join(str(int(bit)) for bit in fractal_output_fixed)
len_fractal_output_fixed = len(fractal_binary_string_fixed)

fractal_binary_string_fixed[:256], len_fractal_output_fixed
