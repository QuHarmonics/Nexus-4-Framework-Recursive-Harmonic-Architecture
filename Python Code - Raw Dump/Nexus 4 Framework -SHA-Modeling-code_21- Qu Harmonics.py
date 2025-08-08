# Natural-form spiral generator based on π digits — no forced torus, uses data-driven curvature

def generate_freeform_pi_spiral(num_digits=10000):
    mp.dps = num_digits + 2
    pi_str = str(mp.pi)[2:num_digits+2]

    # Convert to digits
    pi_digits = [int(d) for d in pi_str]
    diffs = np.diff(pi_digits)
    bit_lengths = [len(bin(d)[2:]) for d in pi_digits[:-1]]

    # Use phase + bit weight to drive curvature
    phi = (1 + np.sqrt(5)) / 2
    theta = np.cumsum([d * (2 * np.pi / phi) for d in pi_digits[:-1]])  # cumulative angle
    r = np.cumsum([bl * 0.05 for bl in bit_lengths])  # cumulative growth by bit-weight

    # Curvature emerges from the recursive interaction
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = np.cumsum([d * 0.05 for d in diffs])  # vertical delta based on tension

    # Plot it as natural flow
    fig = go.Figure(data=go.Scatter3d(
        x=x,
        y=y,
        z=z,
        mode='lines+markers',
        marker=dict(
            size=2,
            color=diffs,
            colorscale='Electric',
            colorbar=dict(title='Δ Tension'),
            opacity=0.85
        ),
        line=dict(color='rgba(120,120,120,0.2)', width=1)
    ))

    fig.update_layout(
        title=f"Emergent Harmonic Spiral of π (First {num_digits} Digits)",
        scene=dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False)),
        height=900,
        width=900,
        showlegend=False
    )

    fig.show()

# Example usage:
generate_freeform_pi_spiral(500000)

