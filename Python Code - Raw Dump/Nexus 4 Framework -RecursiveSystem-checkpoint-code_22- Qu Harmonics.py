import numpy as np
from mpmath import mp
import plotly.graph_objects as go

def stream_pi_digits(n):
    """Generate n digits of pi as a stream of integers."""
    mp.dps = n + 2
    pi_str = str(mp.pi)[2:n+2]
    for ch in pi_str:
        yield int(ch)
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

        
def generate_pi_spiral_streamed(num_digits=100000, 
                                   drive_mode='digit_value', 
                                   chunk_size=10000):
    """Stream-compatible spiral generator for very large digit counts."""
    pi_stream = stream_pi_digits(num_digits)

    x = []
    y = []
    z = []

    prev_digit = next(pi_stream)
    index = 0

    phi = (1 + np.sqrt(5)) / 2

    for current_digit in pi_stream:
        if drive_mode == 'digit_value':
            r = current_digit * 0.1
            theta = index * (2 * np.pi / phi)
        elif drive_mode == 'digit_ratio':
            diff = current_digit - prev_digit
            r = np.abs(diff) * 0.1
            theta = index * (2 * np.pi / phi)
        elif drive_mode == 'odd_even':
            r = 0.5 if current_digit % 2 == 0 else 1.0
            theta = index * (2 * np.pi / phi)
        elif drive_mode == 'binary_len':
            r = len(bin(current_digit)[2:]) * 0.1
            theta = index * (2 * np.pi / phi)

        x.append(r * np.cos(theta))
        y.append(r * np.sin(theta))
        z.append(index * 0.01)  # vertical growth

        prev_digit = current_digit
        index += 1

    fig = go.Figure(data=go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='lines+markers',
        marker=dict(
            size=2,
            color=np.arange(len(x)),
            colorscale='Plasma',
            colorbar=dict(title='Index'),
            opacity=0.8
        ),
        line=dict(color='rgba(120,120,120,0.2)', width=1)
    ))

    fig.update_layout(
        title=f"Emergent Harmonic Spiral of π (First {num_digits} Digits) - {drive_mode.capitalize()} Method",
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)),
        height=900,
        width=900,
        showlegend=False
    )

    fig.show()

# Example usage
generate_pi_spiral_streamed(100000, drive_mode='digit_value')
generate_pi_spiral_streamed(100000, drive_mode='digit_ratio')
generate_pi_spiral_streamed(100000, drive_mode='odd_even')
generate_pi_spiral_streamed(100000, drive_mode='binary_len')