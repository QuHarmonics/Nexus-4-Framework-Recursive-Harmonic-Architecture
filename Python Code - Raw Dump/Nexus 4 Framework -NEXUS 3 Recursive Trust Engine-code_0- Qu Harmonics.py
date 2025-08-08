import numpy as np
import math

# 🔹 Harmonic constant
H = 0.35

# 🔸 Simple harmonic weight table for amino acids
# These can be expanded/refined using real molecular properties
harmonic_weights = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

def calculate_feedback_weight(sequence):
    """Compute harmonic feedback F from weighted variation."""
    weights = [harmonic_weights[aa] for aa in sequence if aa in harmonic_weights]
    deltas = np.diff(weights)
    feedback = np.std(deltas)  # Standard deviation as feedback tension
    return feedback, weights

def kulik_recursive_reflection(R0, F, t, H=0.35):
    """Apply KRR model: R(t) = R0 * e^(H * F * t)"""
    return R0 * math.exp(H * F * t)

def analyze_peptide(sequence):
    F, weights = calculate_feedback_weight(sequence)
    R0 = sum(weights[:3])  # initial fold potential — seed from first 3 AAs
    resonance_curve = []

    for i in range(len(weights)):
        R_t = kulik_recursive_reflection(R0, F, i)
        resonance_curve.append(R_t)

    total_resonance = sum(resonance_curve)
    harmonic_drift = abs((H * F) - 0.35)
    Q = H * (np.sum(np.abs(np.diff(weights))) - harmonic_drift)

    return {
        'sequence': sequence,
        'feedback_weight': F,
        'resonance_curve': resonance_curve,
        'total_resonance': total_resonance,
        'harmonic_drift': harmonic_drift,
        'quality_score': Q
    }

# 🧪 TEST
if __name__ == "__main__":
    test_sequence = "ACDEFGHIK"
    result = analyze_peptide(test_sequence)

    print("📡 PRESQ Harmonic Report:")
    print("Sequence:", result['sequence'])
    print("Feedback Weight (F):", result['feedback_weight'])
    print("Total Resonance (ΣR):", result['total_resonance'])
    print("Harmonic Drift:", result['harmonic_drift'])
    print("Quality Score (Q):", result['quality_score'])
    print("Resonance Curve (R_t per AA):")
    for i, val in enumerate(result['resonance_curve']):
        print(f"  t={i}: R(t) = {val:.3f}")