import numpy as np
import matplotlib.pyplot as plt

# Parameters
time_steps = 1000
states = ["Faces", "Vase", "Superposition"]
state_probabilities = np.array([0.4, 0.4, 0.2])  # Initial probabilities
feedback_strength = 0.1
noise_level = 0.05

# Tracking states over time
state_history = []

# Feedback loop
def feedback_adjustment(current_state, probabilities):
    if current_state == "Faces":
        probabilities[0] += feedback_strength
    elif current_state == "Vase":
        probabilities[1] += feedback_strength
    probabilities[2] -= feedback_strength  # Decrease superposition likelihood
    return probabilities

# Simulation
for t in range(time_steps):
    # Add noise to probabilities
    state_probabilities += np.random.uniform(-noise_level, noise_level, 3)
    state_probabilities = np.clip(state_probabilities, 0, None)  # Ensure no negative probabilities
    state_probabilities /= np.sum(state_probabilities)  # Normalize to sum to 1

    # Choose a state based on probabilities
    current_state = np.random.choice(states, p=state_probabilities)
    state_history.append(current_state)

    # Apply feedback adjustment
    state_probabilities = feedback_adjustment(current_state, state_probabilities)
    state_probabilities = np.clip(state_probabilities, 0, None)  # Ensure no negative probabilities
    state_probabilities /= np.sum(state_probabilities)  # Normalize again after feedback adjustment

# Visualize results
time = np.arange(time_steps)
state_indices = [states.index(state) for state in state_history]

plt.figure(figsize=(12, 6))
plt.plot(time, state_indices, label="State Oscillation")
plt.yticks(ticks=range(len(states)), labels=states)
plt.xlabel("Time Steps")
plt.ylabel("State")
plt.title("Extended Dual-State Oscillation with Feedback")
plt.legend()
plt.show()
