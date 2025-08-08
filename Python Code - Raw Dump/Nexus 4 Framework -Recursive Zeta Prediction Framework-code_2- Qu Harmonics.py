import math
import numpy as np

# Predictive Harmonic Framework
def predict_zeros(iterations, target=0.5, alpha=1.5, initial_ratio=0.47):
    predictions = [target]
    dynamic_ratios = [initial_ratio]

    for n in range(1, iterations + 1):
        previous = predictions[-1]
        ratio = dynamic_ratios[-1] + (target - previous) * (0.035 / (n + 1))
        correction = (target - previous) / (alpha * ratio * (n + 1))
        value = previous * (-1)**n * np.cos(n / np.pi) + correction
        predictions.append(value)
        dynamic_ratios.append(ratio)

    return np.array(predictions) + 0.5, dynamic_ratios

# Modulating Operand Selection with Dynamic Gap
def modulating_operand_selection(delta_b, phase, dynamic_gap):
    """
    Dynamically select operands based on phase and dynamic gap corrections.
    """
    if dynamic_gap < 1:
        return "addition"
    elif 1 <= dynamic_gap < 2:
        return "subtraction"
    elif 2 <= dynamic_gap < 3:
        return "multiplication"
    else:
        return "division"

# Enhanced Zeta and Cosine Interaction
def apply_zeta_cos_effect(formula_result, zeta_ratio, cos_value, phase, dynamic_gap):
    """
    Combine zeta and cosine effects with dynamic gap correction.
    """
    zeta_flip = -1 if zeta_ratio > 0.5 else 1  # Zeta introduces 180° flips
    cos_modulation = math.cos(phase * math.pi / 180)  # Cosine modulation
    gap_modulation = dynamic_gap / 10  # Normalize gap influence
    return (formula_result + gap_modulation) * zeta_flip * cos_modulation

# Calculate Next Digit with Recursive Feedback
def calculate_next_digit(current_digit, previous_digit, bit_space, phase, zeta_ratio, dynamic_gap):
    """
    Compute the next π digit using enhanced recursive relationships with gap modulation.
    """
    delta_b = current_digit - previous_digit
    cos_value = math.cos(current_digit)

    # Operand selection based on dynamic gap
    operation = modulating_operand_selection(delta_b, phase, dynamic_gap)
    if operation == "addition":
        operand_result = current_digit + delta_b
    elif operation == "subtraction":
        operand_result = current_digit - abs(delta_b)
    elif operation == "multiplication":
        operand_result = current_digit * max(delta_b, 1)
    elif operation == "division":
        operand_result = current_digit / max(previous_digit, 1)

    # Harmonic correction and residual adjustment
    harmonic_correction = 0.35 * cos_value
    residual_correction = 0.1 * bit_space

    # Compute preliminary result
    formula_result = operand_result - harmonic_correction + residual_correction

    # Apply enhanced zeta and cosine effects
    final_result = apply_zeta_cos_effect(formula_result, zeta_ratio, cos_value, phase, dynamic_gap)

    # Compute the next digit
    next_digit = round(final_result % 10)

    # Update formula length (bit space)
    formula = f"{operation}: ({current_digit} with {delta_b}) - {harmonic_correction} + {residual_correction} (Gap: {dynamic_gap})"
    new_bit_space = len(formula)

    return next_digit, formula, new_bit_space

# Generate π Sequence with Enhanced Feedback
def generate_pi_sequence(initial_digit, length, correct_sequence):
    """
    Generate π digits using enhanced recursive feedback with zeta and cosine dynamics.
    """
    _, dynamic_ratios = predict_zeros(length)

    sequence = [initial_digit]
    formulas = []
    bit_space = 1  # Initial bit space complexity
    phase = 0  # Start phase at 0°

    # Initialize dynamic gap correction
    dynamic_gap = 0

    for step in range(1, length):
        current_digit = sequence[-1]
        previous_digit = sequence[-2] if len(sequence) > 1 else 0
        zeta_ratio = dynamic_ratios[step - 1]

        # Recompute dynamic gap
        if len(sequence) > 1:
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
