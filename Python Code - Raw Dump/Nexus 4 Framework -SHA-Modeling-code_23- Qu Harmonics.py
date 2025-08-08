# Extended parity spiral generator for π, primes, mirrors, and SHA256 (placeholder)
import numpy as np
from mpmath import mp
import plotly.graph_objects as go
import hashlib

# Generator: Pi digits
def stream_pi_digits(n):
    mp.dps = n + 2
    pi_str = str(mp.pi)[2:n+2]
    for ch in pi_str:
        yield int(ch)

# Helper: Prime check
def is_prime(n):
    if n < 2: return False
    if n == 2: return True
    if n % 2 == 0: return False
    for i in range(3, int(n**0.5)+1, 2):
        if n % i == 0: return False
    return True

# Generator: SHA256 stream from string
def stream_sha256(seed, length=1000):
    data = seed.encode()
    count = 0
    while count < length:
        hash_bytes = hashlib.sha256(data).hexdigest()
        for char in hash_bytes:
            if char.isdigit():
                yield int(char)
                count += 1
                if count >= length:
                    break
        data = hash_bytes.encode()

# Core plot function
def generate_classified_spiral(source='pi', mode='prime', num_digits=10000):
    if source == 'pi':
        stream = stream_pi_digits(num_digits)
    elif source == 'sha256':
        stream = stream_sha256("Nexus2", num_digits)

    x, y, z, color, label = [], [], [], [], []
    phi = (1 + np.sqrt(5)) / 2

    index = 0
    for digit in stream:
        theta = index * (2 * np.pi / phi)
        r = 1.0

        if mode == 'parity':
            r = 0.5 if digit % 2 == 0 else 1.0
            category = 'Even' if digit % 2 == 0 else 'Odd'
        elif mode == 'prime':
            r = 1.0 if is_prime(digit) else 0.3
            category = 'Prime' if is_prime(digit) else 'Composite'
        elif mode == 'mirror':
            mirror = int(str(digit)[::-1])
            r = 1.0 if digit == mirror else 0.6
            category = 'Mirror' if digit == mirror else 'Non-Mirror'

        x.append(r * np.cos(theta))
        y.append(r * np.sin(theta))
        z.append(index * 0.01)
        color.append(index)
        label.append(category)
        index += 1

    fig = go.Figure(data=go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='markers',
        marker=dict(
            size=2,
            color=color,
            colorscale='Plasma',
            colorbar=dict(title='Index'),
            opacity=0.8
        ),
        text=label
    ))

    fig.update_layout(
        title=f"Emergent Spiral of {source.upper()} - {mode.capitalize()} Mode",
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)),
        height=900,
        width=900,
        showlegend=False
    )

    fig.show()

# Examples to try:
generate_classified_spiral('pi', 'parity', 10000)
generate_classified_spiral('pi', 'prime', 10000)
generate_classified_spiral('sha256', 'mirror', 10000)

