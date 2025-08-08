# Re-executing Simulation 1 due to prior connection issue

import hashlib
import numpy as np
import plotly.graph_objs as go
from sympy import symbols, cos, pi, lambdify

# Step 1: Define input symbolic sequence
sequence = "MKQLQARILAVERYLKDQQLLGIWGCSGKLICTTAVPWNAS"  # HIV gp41 segment

# Step 2: SHA-256 collapse function
def sha256_hash(s):
    return hashlib.sha256(s.encode()).hexdigest()

hashed = sha256_hash(sequence)

# Step 3: Convert hash to bitstream and extract Δπ signatures
def hash_to_bit_deltas(hex_hash):
    binary_stream = bin(int(hex_hash, 16))[2:].zfill(256)
    deltas = [abs(int(binary_stream[i]) - int(binary_stream[i+1])) for i in range(len(binary_stream) - 1)]
    return deltas, binary_stream

deltas, binary_stream = hash_to_bit_deltas(hashed)

# Step 4: Define Δπ harmonic model
x = symbols('x')
harmonic_model = cos(2 * pi * x) + cos(4 * pi * x)
harmonic_func = lambdify(x, harmonic_model, 'numpy')

x_vals = np.linspace(0, 1, len(deltas))
y_vals = harmonic_func(x_vals)

# Step 5: Apply resonance comparison
resonance = np.multiply(y_vals, deltas)

# Step 6: Visualization with Plotly
trace1 = go.Scatter(x=x_vals, y=resonance, mode='lines', name='Δπ Resonance Signature')
trace2 = go.Scatter(x=x_vals, y=y_vals, mode='lines', name='Harmonic Baseline', line=dict(dash='dot'))

layout = go.Layout(
    title="Recursive Collapse: Δπ Resonance Model",
    xaxis=dict(title="Normalized Position"),
    yaxis=dict(title="Amplitude"),
    showlegend=True
)

fig = go.Figure(data=[trace1, trace2], layout=layout)
fig.show()

# Return essential data
{
    "original_sequence": sequence,
    "sha256_hash": hashed,
    "resonance_sample": resonance[:10]
}
