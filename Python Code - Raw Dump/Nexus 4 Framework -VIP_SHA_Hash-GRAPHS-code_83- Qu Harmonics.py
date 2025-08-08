import math

def universal_formula(A, B, AX, harmonic_constant=0.35):
    # LEN(C) is defined as the absolute value of the ratio A/B.
    # (If B==0, we return infinity for demonstration purposes.)
    ratio_length = abs(A / B) if B != 0 else float('inf')
    
    # Base term: (A^2 + B^2)
    base_term = A**2 + B**2
    
    # Exponential term: (1 + exp(-10*(AX - harmonic_constant)))
    exp_term = 1 + math.exp(-10 * (AX - harmonic_constant))
    
    # Universal Formula F:
    F = base_term * ratio_length * exp_term
    return F

# Define a set of test values for A, B, and AX.
test_values = [
    (1, 4, 0.35),
    (1, 4, 0.5),
    (2, 3, 0.35),
    (2, 3, 0.5),
    (1, 2, 0.35),
    (1, 2, 0.5),
    (3, 3, 0.35),
    (3, 3, 0.5)
]

print("Universal Formula Test Results:")
for (A, B, AX) in test_values:
    F_val = universal_formula(A, B, AX)
    print(f"A={A}, B={B}, AX={AX} => F={F_val:.6f}")
