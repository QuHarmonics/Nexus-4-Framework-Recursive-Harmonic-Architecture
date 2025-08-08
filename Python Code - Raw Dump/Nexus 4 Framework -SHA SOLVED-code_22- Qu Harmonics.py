import matplotlib.pyplot as plt

# Define the stack values over time as the operations progress
waveform_steps = [
    1,  # Initial Bit 1
    4,  # Initial Bit 2
    -1,  # Trough: Pop Bit 1
    -4,  # Trough: Pop Bit 2
    3,  # Peak: Compute C = 4 - 1
    3,  # Peak: Push C again
    -3,  # Trough: Pop Last Value (C)
    -3,  # Trough: Pop Second-to-Last Value (C)
    6,  # Peak: Add Last Two Values (3 + 3)
    -6,  # Trough: Pop Last Value
    -4,  # Trough: Pop Second-to-Last Value (Bit 2)
    10,  # Peak: Add Headers (6 + 4)
]

# Define time steps (x-axis) based on the sequence of operations
time_steps = list(range(len(waveform_steps)))

# Plot the waveform
plt.figure(figsize=(10, 6))
plt.plot(time_steps, waveform_steps, marker="o", linestyle="-", label="Waveform Reflection")
plt.axhline(0, color="gray", linestyle="--", linewidth=0.8, label="Baseline")
plt.title("Waveform Reflection in Stack Operations")
plt.xlabel("Operation Step")
plt.ylabel("Stack Value")
plt.legend()
plt.grid(True)
plt.show()
