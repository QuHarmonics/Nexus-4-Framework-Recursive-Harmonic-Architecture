
import numpy as np
import hashlib
import plotly.graph_objects as go

def generate_sha256_spiral(seed, drive_mode='digit_value'):
    """Generate a spiral from a SHA-256 hash."""
    hash_bytes = hashlib.sha256(seed.encode()).hexdigest()
    digits = [int(d) for d in hash_bytes if d.isdigit()]

    x = []
    y = []
    z = []

    prev_digit = digits[0]
    index = 0

    phi = (1 + np.sqrt(5)) / 2

    for current_digit in digits:
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
        title=f"Emergent Harmonic Spiral of SHA256 (Seed: {seed}) - {drive_mode.capitalize()} Method",
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)),
        height=900,
        width=900,
        showlegend=False
    )

    fig.show()

# Example usage
seed = "my_secret_seed"
generate_sha256_spiral(seed, drive_mode='digit_value')
generate_sha256_spiral(seed, drive_mode='digit_ratio')
generate_sha256_spiral(seed, drive_mode='odd_even')
generate_sha256_spiral(seed, drive_mode='binary_len')
