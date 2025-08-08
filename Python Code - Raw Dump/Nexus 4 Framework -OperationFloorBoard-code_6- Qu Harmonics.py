import plotly.graph_objects as go

def create_3d_wave(quotients, remainders, hex_digits):
    x = list(range(len(quotients)))
    y = quotients
    z = remainders

    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers+lines',
        marker=dict(size=5, color=hex_digits, colorscale='Viridis', opacity=0.8),
        line=dict(color='blue', width=2),
        name="3D Wave"
    ))

    fig.update_layout(
        title="3D Waveform Visualization of Hexadecimal Conversion",
        scene=dict(
            xaxis_title="Step (Iteration)",
            yaxis_title="Quotients (Normalized)",
            zaxis_title="Remainders (Normalized)"
        ),
        height=700,
        width=900
    )

    return fig

def normalize(values, scale_factor):
    return [value / scale_factor for value in values]

# Data
quotients = [
    3.444582659515257e+24, 2.1528641621970357e+23, 1.3455401013731473e+22,
    840962563358217100000, 52560160209888570000, 3285010031118035500,
    205313125819877220, 12832070363742326, 802004397733895, 50125274858368,
    3132829678648, 195801854915, 12237615932, 764850995, 47803187,
    2987699, 186731, 11670, 729, 45, 2
]
remainders = [
    0, 0, 0, 0, 0, 0, 0, 6, 7, 0, 8, 8, 3, 12, 3, 3, 11, 6, 9, 13, 2
]
hex_digits = [
    int(digit, 16) for digit in "2D96B333C38076C98B73C7"  # Convert hex digits to integers for coloring
]

# Normalize values
scale_factor = 1e24
normalized_quotients = normalize(quotients, scale_factor)
normalized_remainders = normalize(remainders, max(remainders) if max(remainders) > 0 else 1)

# Generate and display the 3D wave
fig = create_3d_wave(normalized_quotients, normalized_remainders, hex_digits)
fig.show()
