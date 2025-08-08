import numpy as np

def projective_coordinate(D):
    """Normalize vector D in Z9^8 to P7(Z9)."""
    g = np.gcd.reduce(D)
    return (D // g) % 9

def classify(D, lut):
    coord = tuple(projective_coordinate(D))
    return lut[coord]
