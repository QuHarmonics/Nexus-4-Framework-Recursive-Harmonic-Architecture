import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

def strain_to_hex(s):
    # s in [0,1]
    if s <= 0.5:
        t = 2*s
        R = (1-t)*0 + t*170
        G = (1-t)*0 + t*170
        B = (1-t)*255 + t*170
    else:
        t = 2*s - 1
        R = (1-t)*170 + t*255
        G = (1-t)*170 + t*0
        B = (1-t)*170 + t*0
    return f"#{int(R):02X}{int(G):02X}{int(B):02X}"

# Example hexagon lattice coords (hex_spiral from before)
# and sample strain values epsilon_i
coords = hex_spiral(64)  # axial coords
xs = [np.sqrt(3)*(q + r/2) for q,r in coords]
ys = [1.5*r for q,r in coords]

# Dummy strain field for demonstration (e.g. bending around x=0)
xs_arr = np.array(xs)
eps = (xs_arr - xs_arr.min())/(xs_arr.max()-xs_arr.min())  # simple gradient
strain_vals = eps  # normalized 0–1

fig, ax = plt.subplots(figsize=(6,6))
R = 1.0
for x, y, s in zip(xs, ys, strain_vals):
    color = strain_to_hex(s)
    hexagon = RegularPolygon((x, y), numVertices=6, radius=R,
                             orientation=np.pi/6, facecolor=color, edgecolor='gray')
    ax.add_patch(hexagon)
ax.set_aspect('equal')
ax.axis('off')
plt.show()