# ─────────────────────────────────────────────────────────────────────────────
#  Quantum Pi Estimator — Optimized Full Version
# ─────────────────────────────────────────────────────────────────────────────

import numpy as np
import matplotlib.pyplot as plt
from qiskit import QuantumCircuit, Aer, execute
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import QFT

# Create Quantum Phase Estimation Circuit
def run_quantum_phase_estimation(bits, shots=4096):
    qc = QuantumCircuit(bits + 1, bits)

    # Prepare the target qubit in |1> state
    qc.x(bits)

    # Apply Hadamard gates to the counting qubits
    for qubit in range(bits):
        qc.h(qubit)

    # Controlled phase shifts (simulate U with eigenvalue phase = pi/4)
    repetitions = 1
    for counting_qubit in range(bits):
        for _ in range(repetitions):
            qc.cp(np.pi, counting_qubit, bits)
        repetitions *= 2

    # Apply inverse QFT
    qc.append(QFT(num_qubits=bits, do_swaps=True).inverse(), range(bits))

    # Measurement
    qc.measure(range(bits), range(bits))

    # Execute on Aer simulator
    simulator = Aer.get_backend('qasm_simulator')
    job = execute(qc, simulator, shots=shots)
    result = job.result()
    counts = result.get_counts(qc)
    
    return counts

# Refine Pi estimate based on measurements
def refine_pi_estimate(previous_estimate, new_measurement, learning_rate=0.2):
    measured_phase = int(new_measurement, 2) / (2 ** len(new_measurement))
    updated_pi = previous_estimate + learning_rate * ((measured_phase * 4) - previous_estimate)
    return updated_pi

# Initialize
initial_pi_estimate = 3.0  # Starting point
pi_estimates = [initial_pi_estimate]
bit_precision = 5          # Precision bits
iterations = 20            # How many times to refine

# Perform the iterative estimation
for _ in range(iterations):
    counts = run_quantum_phase_estimation(bit_precision)
    most_probable_result = max(counts, key=counts.get)  # Get most frequent result
    updated_pi = refine_pi_estimate(pi_estimates[-1], most_probable_result)
    pi_estimates.append(updated_pi)

# Plot Convergence
plt.figure(figsize=(10, 6))
plt.plot(range(len(pi_estimates)), pi_estimates, marker='o', linestyle='-', color='blue')
plt.axhline(y=np.pi, color='red', linestyle='--', label='True Pi ≈ 3.14159...')
plt.title('Quantum Phase Estimation: Convergence to Pi')
plt.xlabel('Iteration')
plt.ylabel('Estimate of Pi')
plt.legend()
plt.grid(True)
plt.show()

# Display final estimate
print(f"\nFinal estimated π value after {iterations} iterations: {pi_estimates[-1]}")
print(f"Real π value: {np.pi}")
print(f"Absolute error: {abs(pi_estimates[-1] - np.pi)}")
