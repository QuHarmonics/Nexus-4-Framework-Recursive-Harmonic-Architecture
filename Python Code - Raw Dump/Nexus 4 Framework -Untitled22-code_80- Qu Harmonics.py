import numpy as np
from scipy.special import zeta

def elliptic_curve_discriminant(a, b):
    return -16 * (4 * a**3 + 27 * b**2)

def l_function_numerical(a, b, s=1.0):
    discriminant = elliptic_curve_discriminant(a, b)
    if discriminant == 0:
        return None  # Skip calculation if the curve is singular
    return zeta(s) * np.exp(-np.abs(a + b) * s)  # Simplified L-function model

def simulate_elliptic_curves(a_range, b_range, num_points):
    results = []
    for a in np.linspace(*a_range, num=num_points):
        for b in np.linspace(*b_range, num=num_points):
            if elliptic_curve_discriminant(a, b) != 0:
                l_value = l_function_numerical(a, b)
                if l_value is not None:
                    results.append({'a': a, 'b': b, 'L(s=1)': l_value})
    return results

# Adjusted parameter ranges to potentially reduce the chance of singularity
a_range = (-3, 3)  # Increased range for broader curve generation
b_range = (-3, 3)
num_points = 50  # Number of points in each range

# Run the simulation with the adjusted ranges
simulation_results = simulate_elliptic_curves(a_range, b_range, num_points)
for result in simulation_results:
    print(result)
