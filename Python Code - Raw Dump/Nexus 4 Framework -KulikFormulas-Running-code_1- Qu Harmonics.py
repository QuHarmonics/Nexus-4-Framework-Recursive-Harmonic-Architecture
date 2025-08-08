import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque
from mpmath import mp

# ------------------- Precision & π Configuration ------------------- #
mp.dps = 10000  # Set high-precision digits for π
pi_digits = str(mp.pi)[2:]  # Skip "3."

# ------------------- Layer 1: π Digit Stream as Entropic Feed ------------------- #
def get_pi_window(offset, length=8):
    digits = pi_digits[offset:offset + length]
    return [int(d) for d in digits if d.isdigit()]

# ------------------- Layer 2: Harmonic Feedback Growth ------------------- #
def harmonic_growth_step(seq, target_h=0.35):
    h = sum(seq) / len(seq)
    deviation = abs(h - target_h)
    if deviation < 0.01:
        return None, h
    next_val = (seq[-1] + seq[-2]) % 10
    return next_val, h

# ------------------- Live Memory Buffer ------------------- #
window = 120  # Display length
seq = deque([3, 2], maxlen=window)
h_values = deque(maxlen=window)
x_vals = deque(range(window), maxlen=window)
pi_offset = 0

# ------------------- Plot Setup ------------------- #
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6))
line1, = ax1.plot([], [], lw=2, color='orange', label='Heartbeat')
line2, = ax2.plot([], [], lw=1.5, color='blue', label='H(t)')
ax1.set_ylim(0, 10)
ax2.set_ylim(0, 2)
ax1.set_xlim(0, window)
ax2.set_xlim(0, window)
ax1.set_title('Real-Time Recursive Harmonic Heartbeat')
ax1.set_ylabel('Byte Value')
ax2.set_ylabel('H Value')
ax2.axhline(0.35, color='gray', linestyle='--', lw=0.8)
ax2.set_xlabel('Time')
ax1.grid(True)
ax2.grid(True)
ax1.legend()
ax2.legend()

# ------------------- Update Logic ------------------- #
def update(frame):
    global pi_offset, seq, h_values

    next_val, h = harmonic_growth_step(list(seq))
    if next_val is None or (seq[-1] == 0 and seq[-2] == 0):
        # System collapse: inject fresh entropy using next π window
        pi_seed = get_pi_window(pi_offset, length=8)
        pi_offset = (pi_offset + 7) % len(pi_digits)
        if len(pi_seed) < 2:
            pi_seed = [3, 2]
        seq.extend(pi_seed)
    else:
        seq.append(next_val)

    h_values.append(h if h is not None else 0)
    x_vals.append(x_vals[-1] + 1)

    # Update plots
    line1.set_data(range(len(seq)), list(seq))
    line2.set_data(range(len(h_values)), list(h_values))

    ax1.set_xlim(max(0, x_vals[-1] - window), x_vals[-1])
    ax2.set_xlim(max(0, x_vals[-1] - window), x_vals[-1])
    return line1, line2

# ------------------- Start Animation Loop ------------------- #
ani = animation.FuncAnimation(
    fig,
    update,
    interval=200,
    blit=True,
    cache_frame_data=False  # ✅ Avoid unbounded memory growth
)
plt.tight_layout()
plt.show()
