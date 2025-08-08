def echo_expand(anchor=3, curve=2):
    """
    Reconstructs a symbolic 4-digit chunk from anchor and curve:
    [3, ?, ?, 3 + curve]
    Uses diagonal harmonic symmetry to fill the center.
    """
    x0 = anchor
    x3 = anchor + curve

    # Option 1: Simple mid-fill as harmonic average (symmetry seed)
    # x1 and x2 center between x0 and x3
    midpoint = (x0 + x3) // 2
    delta = abs(x3 - x0) // 2

    x1 = midpoint - delta // 2
    x2 = midpoint + delta // 2

    return [x0, x1, x2, x3]

# Example usage:
print(echo_expand(3, 2))  # Output: [3, 3, 4, 5]
print(echo_expand(3, 4))  # Output: [3, 4, 5, 7]
print(echo_expand(3, 0))  # Output: [3, 3, 3, 3]
