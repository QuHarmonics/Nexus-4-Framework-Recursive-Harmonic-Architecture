import matplotlib.animation as animation

def update_path(frame, path, plot):
    plot.set_data(path[:frame, 0], path[:frame, 1])
    return plot,

fig, ax = plt.subplots()
path = np.array([[np.cos(t), np.sin(t)] for t in np.linspace(0, 0.35*2*np.pi, 100)])
plot, = ax.plot([], [], 'r-')
ani = animation.FuncAnimation(fig, update_path, frames=100, fargs=(path, plot))
plt.show()