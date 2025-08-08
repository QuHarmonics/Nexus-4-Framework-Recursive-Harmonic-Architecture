import math

A = 3
B = 4
C_len = 2
Ax = 0.5

# Calculate each part of the formula
geometric_part = A**2 + B**2
dynamic_capacity = C_len
nonlinear_adjustment = 1 + math.exp(-10 * (Ax - 0.35))

# Final computation
F = geometric_part * dynamic_capacity * nonlinear_adjustment
F
