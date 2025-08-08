import matplotlib.pyplot as plt
import numpy as np

hex_digits = [2, 4, 3, 15, 6, 10, 8, 8, 8, 5]
angles = [np.pi / 8 * d for d in hex_digits]
x, y = 0, 0
xs, ys = [x], [y]

for angle in angles:
    x += np.cos(angle)
    y += np.sin(angle)
    xs.append(x)
    ys.append(y)

plt.plot(xs, ys, marker='o')
plt.axis('equal')
plt.title('BBP Hex Digits as Angular Walk')
plt.show()
