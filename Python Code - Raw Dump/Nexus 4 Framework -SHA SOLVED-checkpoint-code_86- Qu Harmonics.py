import numpy as np
import matplotlib.pyplot as plt
import random
import time

def simulate_dual_wave_states(total_steps=200, observation_interval=20):
    """
    Simulates a dual-state wave system oscillating between two interpretations (e.g., faces and vase).

    Parameters:
    - total_steps: Total number of simulation steps.
    - observation_interval: Number of steps before a blink (observation event).

    Returns:
    - states: The states over time ("Faces" or "Vase").
    - collapsed_states: The states when observed (collapsed state).
    """
    states = []
    collapsed_states = []

    # Dual-state probabilities
    state_probabilities = {"Faces": 0.5, "Vase": 0.5}
    current_state = random.choice(["Faces", "Vase"])

    for step in range(total_steps):
        # Oscillation: Randomly switch between states based on probabilities
        if random.random() < state_probabilities["Faces"]:
            current_state = "Faces"
        else:
            current_state = "Vase"
        states.append(current_state)

        # Observation (blink) forces collapse to one state
        if step % observation_interval == 0:
            collapsed_state = random.choice(["Faces", "Vase"])
            collapsed_states.append((step, collapsed_state))

    return states, collapsed_states

def visualize_dual_wave_states(states, collapsed_states, total_steps=200):
    """
    Visualizes the dual-state wave system and its observed (collapsed) states.

    Parameters:
    - states: List of states over time.
    - collapsed_states: List of collapsed states due to observation.
    - total_steps: Total number of simulation steps.
    """
    time_steps = np.arange(total_steps)
    state_mapping = {"Faces": 1, "Vase": -1}
    state_values = [state_mapping[state] for state in states]

    plt.figure(figsize=(10, 6))

    # Plot the oscillating dual states
    plt.plot(time_steps, state_values, label="Oscillating States", alpha=0.7)

    # Highlight the collapsed states
    for step, state in collapsed_states:
        plt.axvline(x=step, color='red', linestyle='--', alpha=0.5)
        plt.text(step, state_mapping[state], state, color='red', fontsize=10)

    plt.title("Simulation of Dual-State Oscillation and Collapse")
    plt.xlabel("Time Steps")
    plt.ylabel("State (-1 = Vase, 1 = Faces)")
    plt.yticks([-1, 1], ["Vase", "Faces"])
    plt.legend()
    plt.grid(True)
    plt.show()

def main():
    total_steps = 200
    observation_interval = 20

    print("Simulating dual-state oscillation and collapse...")
    states, collapsed_states = simulate_dual_wave_states(total_steps, observation_interval)
    print(f"Simulation complete. Total steps: {total_steps}, Observation interval: {observation_interval}")

    visualize_dual_wave_states(states, collapsed_states, total_steps)

if __name__ == "__main__":
    main()
