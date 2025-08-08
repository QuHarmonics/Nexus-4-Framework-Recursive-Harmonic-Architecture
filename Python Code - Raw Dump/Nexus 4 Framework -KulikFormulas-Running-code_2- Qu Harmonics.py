import matplotlib.pyplot as plt
import matplotlib.animation as animation
from collections import deque
from mpmath import mp

# Set high-precision π
mp.dps = 10000
pi_digits = str(mp.pi)[2:]

# --- π digit extractor ---
def get_pi_window(offset, length=8):
    digits = pi_digits[offset:offset + length]
    return [int(d) for d in digits if d.isdigit()][:length]

# --- Harmonic feedback step ---
def harmonic_growth_step(seq, target_h=0.35):
    h = sum(seq) / len(seq)
    deviation = abs(h - target_h)
    if deviation < 0.01:
        return None, h
    next_val = (seq[-1] + seq[-2]) % 10
    return next_val, h

# --- Sliding window config ---
window = 200
seq = deque([3, 2], maxlen=window)
h_vals = deque([0.0, 0.0], maxlen=window)
pi_offset = 0
frame_idx = 0

# --- Setup plot ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))
line1, = ax1.plot([], [], lw=2, color='orange', label='Heartbeat')
line2, = ax2.plot([], [], lw=1.5, color='blue', label='H(t)')

# Configure plot ranges
ax1.set_ylim(0, 10)
ax2.set_ylim(0, 2)
ax2.axhline(0.35, color='gray', linestyle='--', lw=0.8)
ax1.set_ylabel('Byte Value')
ax2.set_ylabel('H Value')
ax2.set_xlabel('Time')
ax1.grid(True)
ax2.grid(True)
ax1.legend()
ax2.legend()

# --- Update function ---
def update(frame):
    global pi_offset, seq, h_vals, frame_idx
    frame_idx += 1

    next_val, h = harmonic_growth_step(list(seq))
    if next_val is None or (seq[-1] == 0 and seq[-2] == 0):
        # Collapse detected: reseed with fresh π digits
        pi_seed = get_pi_window(pi_offset, length=8)
        pi_offset = (pi_offset + 7) % len(pi_digits)
        if len(pi_seed) >= 2:
            for d in pi_seed:
                seq.append(d)
                h_vals.append(sum(seq) / len(seq))
    else:
        seq.append(next_val)
        h_vals.append(h if h is not None else 0)

    x = list(range(frame_idx - len(seq) + 1, frame_idx + 1))

    # Update plot data
    line1.set_data(x, list(seq))
    line2.set_data(x, list(h_vals))
    ax1.set_xlim(x[0], x[-1])
    ax2.set_xlim(x[0], x[-1])

    return line1, line2

# --- Animate ---
ani = animation.FuncAnimation(
    fig, update, interval=200, blit=True, cache_frame_data=False
)

plt.show()
