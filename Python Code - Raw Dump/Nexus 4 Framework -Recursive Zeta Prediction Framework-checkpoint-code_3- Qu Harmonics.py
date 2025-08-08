def generate_pi_sequence(initial_digit, length, correct_sequence):
    """
    Generate π digits using enhanced recursive feedback with zeta and cosine dynamics.
    """
    _, dynamic_ratios = predict_zeros(length)

    sequence = [initial_digit]
    formulas = []
    bit_space = 1  # Initial bit space complexity
    phase = 0  # Start phase at 0°

    # Extend correct_sequence if it's shorter than the desired length
    if len(correct_sequence) < length:
        correct_sequence = correct_sequence + [0] * (length - len(correct_sequence))

    # Initialize dynamic gap correction
    dynamic_gap = 0

    for step in range(1, length):
        current_digit = sequence[-1]
        previous_digit = sequence[-2] if len(sequence) > 1 else 0
        zeta_ratio = dynamic_ratios[step - 1]

        # Recompute dynamic gap
        dynamic_gap = abs(correct_sequence[step - 1] - sequence[-1])

        next_digit, formula, new_bit_space = calculate_next_digit(
            current_digit, previous_digit, bit_space, phase, zeta_ratio, dynamic_gap
        )
        sequence.append(next_digit)
        formulas.append(formula)

        # Update bit space and phase
        bit_space = new_bit_space
        phase = (phase + 90) % 360  # Increment phase in 90° steps

    return sequence, formulas

# Example Usage
initial_digit = 3
sequence_length = 50
correct_sequence = [3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 8, 9, 7, 9, 3, 2, 3, 8, 4]  # Correct π sequence snippet
pi_sequence, formulas_used = generate_pi_sequence(initial_digit, sequence_length, correct_sequence)

print("Generated π Sequence:", pi_sequence)
print("\nFormulas Used:")
for idx, formula in enumerate(formulas_used, start=1):
    print(f"Step {idx}: {formula}")
