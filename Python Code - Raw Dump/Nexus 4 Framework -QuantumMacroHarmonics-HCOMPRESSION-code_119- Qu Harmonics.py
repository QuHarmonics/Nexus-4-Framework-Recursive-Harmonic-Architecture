import numpy as np

def linear_feedback_control(x, desired_state, feedback_gain, alpha=0.1):
    harmonic_term = alpha * x * np.cos(x)
    return feedback_gain * (desired_state - x) + harmonic_term

def nonlinear_feedback_control(x, desired_state, feedback_gain, alpha=0.1):
    harmonic_term = alpha * x * np.cos(x)
    return feedback_gain * (desired_state - x) ** 3 + harmonic_term

# Lorenz-84 system parameters
x0, y0, z0 = (0, 1, 1)

# Linear feedback control with harmonic feedback
x_linear_harmonic = x0
tolerance = 1e-6
while True:
    x_new = x_linear_harmonic + linear_feedback_control(x_linear_harmonic, 0, 0.1)
    if abs(x_new - x_linear_harmonic) < tolerance:
        break
    x_linear_harmonic = x_new

# Nonlinear feedback control with harmonic feedback
x_nonlinear_harmonic = x0
tolerance = 1e-6
while True:
    x_new = x_nonlinear_harmonic + nonlinear_feedback_control(x_nonlinear_harmonic, 0, 0.1)
    if abs(x_new - x_nonlinear_harmonic) < tolerance:
        break
    x_nonlinear_harmonic = x_new

print("Linear feedback control with harmonic feedback:", x_linear_harmonic)
print("Nonlinear feedback control with harmonic feedback:", x_nonlinear_harmonic)