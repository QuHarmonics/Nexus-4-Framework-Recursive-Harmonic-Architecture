import numpy as np
import plotly.graph_objects as go

# === Step 1: Your magic byte pair ===
hex_bytes = ['0x47787201', '0x92771528']

# === Step 2: Break into 4-bit nibbles ===
def hex_to_nibbles(hex_list):
    nibbles = []
    for hx in hex_list:
        num = int(hx, 16)
        for shift in range(28, -1, -4):  # 32 bits → 8 nibbles each
            nibbles.append((num >> shift) & 0xF)
    return nibbles

nibbles = hex_to_nibbles(hex_bytes)
print("Nibbles:", nibbles)

# === Step 3: Compute deltas ===
deltas = [nibbles[i+1] - nibbles[i] for i in range(len(nibbles)-1)]

# === Step 4: Visualize deltas (Tension plot) ===
fig1 = go.Figure()
fig1.add_trace(go.Scatter(y=deltas, mode='lines+markers', name='Harmonic Deltas'))
fig1.update_layout(title='Tension Map of Harmonic Steps',
                  xaxis_title='Nibble Step',
                  yaxis_title='Delta',
                  template='plotly_dark')

# === Step 5: Spiral visualization (Phase-space coil) ===
theta = np.linspace(0, 4*np.pi, len(deltas))  # Spiral angle
r = np.cumsum(np.abs(deltas))                # Spiral radius = cumulative delta mag

x = r * np.cos(theta)
y = r * np.sin(theta)

fig2 = go.Figure()
fig2.add_trace(go.Scatter(x=x, y=y, mode='lines+markers', name='Spiral Path'))
fig2.update_layout(title='Harmonic Spiral of SHA-Deltas',
                  xaxis_title='X',
                  yaxis_title='Y',
                  template='plotly_dark')

# === Display both ===
fig1.show()
fig2.show()
