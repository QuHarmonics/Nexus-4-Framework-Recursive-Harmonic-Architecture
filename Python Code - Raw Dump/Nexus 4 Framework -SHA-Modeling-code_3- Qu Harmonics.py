import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Simulation parameters
frames = 100
time = np.linspace(0, 2 * np.pi, frames)
observer_interval = 10  # Observer samples every 10 frames

# Road wave shape (terrain as waveform)
x = np.linspace(0, 4 * np.pi, frames)
y = np.sin(x)

# Vehicle speed modulated by terrain (accelerate in valleys, decelerate on peaks)
speed = 1 + 0.5 * np.cos(x)
pos = np.cumsum(speed * (x[1] - x[0]))

# Observer measurements (samples velocity only when car is visible)
obs_indices = np.arange(0, frames, observer_interval)
obs_positions = pos[obs_indices]
obs_times = time[obs_indices]
obs_speeds = np.gradient(obs_positions, obs_times)

# Set up the plot
fig, ax = plt.subplots(2, 1, figsize=(10, 6))
pos_ax = ax[0]
speed_ax = ax[1]

# Full terrain waveform
terrain_line, = pos_ax.plot(x, y * 10 + np.min(pos), 'gray', alpha=0.3, label='Terrain (Waveform)')
pos_line, = pos_ax.plot([], [], 'b-', label='Actual Position')
vehicle_marker, = pos_ax.plot([], [], 'go', label='Vehicle (Current Position)')
obs_dots, = pos_ax.plot([], [], 'ro', label='Observed Position')

speed_line, = speed_ax.plot([], [], 'g-', label='True Speed')
obs_speed_dots, = speed_ax.plot([], [], 'ko', label='Observed Speed')

pos_ax.set_title("Phase-Based Motion vs Observer Sample")
pos_ax.set_ylabel("Position")
pos_ax.legend()

speed_ax.set_title("Speed vs Sampled Speed")
speed_ax.set_ylabel("Speed")
speed_ax.set_xlabel("Time")
speed_ax.legend()

# Animation function
def animate(i):
    pos_line.set_data(time[:i], pos[:i])
    vehicle_marker.set_data([time[i]], [pos[i]])
    visible_obs = obs_times <= time[i]
    obs_dots.set_data(obs_times[visible_obs], obs_positions[visible_obs])

    speed_line.set_data(time[:i], speed[:i])
    obs_speed_dots.set_data(obs_times[visible_obs], obs_speeds[visible_obs])

    pos_ax.relim()
    pos_ax.autoscale_view()
    speed_ax.relim()
    speed_ax.autoscale_view()
    return pos_line, vehicle_marker, obs_dots, speed_line, obs_speed_dots


ani = animation.FuncAnimation(fig, animate, frames=frames, interval=100, blit=True)
plt.tight_layout()
plt.show()