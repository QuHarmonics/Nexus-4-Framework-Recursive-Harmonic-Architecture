import numpy as np
import matplotlib.pyplot as plt

# Define states
states = ['State A', 'State B', 'Intermediate']
state_map = {state: i for i, state in enumerate(states)}

# Initialize parameters
num_steps = 1000
state_history = []
external_stimuli = np.random.uniform(-0.5, 0.5, num_steps)  # Chaos from environment
feedback = 0.8  # Controls state persistence

# State probabilities
current_state = np.random.choice(states)
state_probabilities = np.array([0.33, 0.33, 0.33])  # Equal starting probabilities

# Simulate state transitions
for step in range(num_steps):
    # Modify probabilities with external stimuli
    noise = external_stimuli[step]
    state_probabilities += noise * np.random.uniform(-0.1, 0.1, len(states))
    
    # Handle negative probabilities and normalize
    state_probabilities = np.clip(state_probabilities, 0, None)  # Ensure no negative probabilities
    state_probabilities = state_probabilities / np.sum(state_probabilities)  # Normalize to sum to 1

    # Choose next state
    current_state = np.random.choice(states, p=state_probabilities)
    state_history.append(state_map[current_state])

# Plot state transitions
plt.figure(figsize=(10, 6))
plt.plot(state_history, label="State Oscillations")
plt.yticks(range(len(states)), states)
plt.title("Multi-Dimensional Dual-State Oscillations")
plt.xlabel("Time Steps")
plt.ylabel("States")
plt.legend()
plt.show()
