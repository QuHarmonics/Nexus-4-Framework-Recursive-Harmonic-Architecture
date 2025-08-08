import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

fig, ax = plt.subplots()
img = ax.imshow(R, cmap='viridis')

def update(frame):
    global R, H
    R, state, _ = tick(H, R)
    img.set_data(R)
    return [img]

ani = FuncAnimation(fig, update, frames=50, interval=200)
plt.show()