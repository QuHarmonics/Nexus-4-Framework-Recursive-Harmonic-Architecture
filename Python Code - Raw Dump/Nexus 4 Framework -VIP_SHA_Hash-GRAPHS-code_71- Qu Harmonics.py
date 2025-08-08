###############################################################################
# Quantum Pi Estimation — FULL ANIMATED VERSION
###############################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from qiskit import QuantumCircuit, Aer, execute
from qiskit.circuit.library import QFT

# ─────────────────────────────────────────────────────────────────────────────
# 1. Core Quantum Phase Estimation functions
# ─────────────────────────────────────────────────────────────────────────────

def run_quantum_phase_estimation(bits, shots=4096):
    qc = QuantumCircuit(bits + 1, bits)

    # Prepare the qubit in |1> state
    qc.x(bits)

    # Apply Hadamard gates
    for qubit in range(bits):
        qc.h(qubit)

    # Controlled phase rotations
    repetitions = 1
    for counting_qubit in range(bits):
        for _ in range(repetitions):
            qc.cp(np.pi, counting_qubit, bits)
        repetitions *= 2

    # Inverse QFT
    qc.append(QFT(num_qubits=bits, do_swaps=True).inverse(), range(bits))

    # Measure
    qc.measure(range(bits), range(bits))

    simulator = Aer.get_backend('qasm_simulator')
    job = execute(qc, simulator, shots=shots)
    result = job.result()
    counts = result.get_counts(qc)

    return counts

def refine_pi_estimate(previous_estimate, new_measurement, learning_rate=0.2):
    measured_phase = int(new_measurement, 2) / (2 ** len(new_measurement))
    updated_pi = previous_estimate + learning_rate * ((measured_phase * 4) - previous_estimate)
    return updated_pi

# ─────────────────────────────────────────────────────────────────────────────
# 2. Animation Setup
# ─────────────────────────────────────────────────────────────────────────────

# Initial parameters
initial_pi_estimate = 3.0
pi_estimates = [initial_pi_estimate]
bit_precision = 5
iterations = 20

# Storage
frames = []

# Pre-run the estimation to store frames
current_estimate = initial_pi_estimate
for _ in range(iterations):
    counts = run_quantum_phase_estimation(bit_precision)
    most_probable_result = max(counts, key=counts.get)
    current_estimate = refine_pi_estimate(current_estimate, most_probable_result)
    pi_estimates.append(current_estimate)
    frames.append(list(pi_estimates))  # store a copy of current history

# ─────────────────────────────────────────────────────────────────────────────
# 3. Animate the convergence
# ─────────────────────────────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(10, 6))
true_pi_line, = ax.plot([], [], 'r--', label='True π ≈ 3.14159...')
estimate_line, = ax.plot([], [], 'bo-', label='Estimated π')
ax.set_xlim(0, iterations)
ax.set_ylim(2.5, 3.5)
ax.set_title('Quantum Phase Estimation: Animated Convergence to Pi')
ax.set_xlabel('Iteration')
ax.set_ylabel('Estimate of Pi')
ax.grid(True)
ax.legend()

def init():
    estimate_line.set_data([], [])
    return estimate_line, true_pi_line

def update(frame_idx):
    current_x = list(range(len(frames[frame_idx])))
    current_y = frames[frame_idx]
    estimate_line.set_data(current_x, current_y)
    return estimate_line, true_pi_line

ani = FuncAnimation(fig, update, frames=len(frames), init_func=init,
                    blit=True, interval=500, repeat=True)

plt.show()

# ─────────────────────────────────────────────────────────────────────────────
# 4. Final Output Summary
# ─────────────────────────────────────────────────────────────────────────────

final_estimate = pi_estimates[-1]
print(f"\nFinal estimated π after {iterations} iterations: {final_estimate}")
print(f"Real π value: {np.pi}")
print(f"Absolute error: {abs(final_estimate - np.pi)}")
