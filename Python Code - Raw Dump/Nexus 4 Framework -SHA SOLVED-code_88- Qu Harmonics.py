import matplotlib.pyplot as plt
import numpy as np

# Parameters
time_steps = 500
observation_interval = 50
noise_level = 0.1
states = ["Faces", "Vase", "Superposition"]
probabilities = [0.5, 0.5, 0.0]  # Initial probabilities

# Initialize state history
state_history = []

# Function to collapse state
def collapse_state(probs):
    max_prob = max(probs)
    return np.random.choice(states, p=[p / sum(probs) for p in probs])

# Simulate dual-state oscillation with noise
for t in range(time_steps):
    # Add noise to probabilities
    noisy_probs = [p + np.random.uniform(-noise_level, noise_level) for p in probabilities]
    noisy_probs = [max(0, p) for p in noisy_probs]  # Ensure no negative probabilities
    noisy_probs = [p / sum(noisy_probs) for p in noisy_probs]  # Normalize

    # Oscillate probabilities
    probabilities[0] = 0.5 + 0.5 * np.sin(2 * np.pi * t / observation_interval)
    probabilities[1] = 1 - probabilities[0]

    if t % observation_interval == 0:
        # Collapse state at observation intervals
        collapsed_state = collapse_state(noisy_probs)
        probabilities = [0.0, 0.0, 0.0]  # Reset probabilities after collapse
        probabilities[states.index(collapsed_state)] = 1.0
        state_history.append((t, collapsed_state))
    else:
        state_history.append((t, "Superposition"))

# Extract time and states
time, observed_states = zip(*state_history)

# Plot results
plt.figure(figsize=(12, 6))
plt.title("Extended Dual-State Oscillation and Collapse")
plt.plot(time, [states.index(s) for s in observed_states], label="State Oscillation", color="blue", alpha=0.6)
plt.xlabel("Time Steps")
plt.ylabel("State")
plt.yticks(range(len(states)), states)
plt.legend()
plt.grid()
plt.show()
