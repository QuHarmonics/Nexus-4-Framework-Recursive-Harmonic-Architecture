from qiskit import QuantumCircuit, Aer, execute, transpile, assemble
from qiskit.visualization import plot_histogram
from qiskit.circuit.library import QFT
from math import pi
import matplotlib.pyplot as plt

def run_quantum_phase_estimation(bits, shots=1024):
    # Create a Quantum Circuit with enough qubits
    qc = QuantumCircuit(bits + 1, bits)
    
    # Prepare the qubit in the |1> state
    qc.x(bits)

    # Apply Hadamard gates to the first n qubits
    for qubit in range(bits):
        qc.h(qubit)
    
    # Controlled-U operations
    repetitions = 1
    for counting_qubit in range(bits):
        for _ in range(repetitions):
            qc.cp(pi/2**counting_qubit, counting_qubit, bits)  # Controlled-phase gate
        repetitions *= 2

    # Apply inverse QFT
    qc.append(QFT(num_qubits=bits, do_swaps=False).inverse(), range(bits))

    # Measure
    for n in range(bits):
        qc.measure(n, n)
    
    # Execute the circuit on the qasm simulator
    simulator = Aer.get_backend('qasm_simulator')
    transpiled_qc = transpile(qc, simulator)
    qobj = assemble(transpiled_qc, shots=shots)
    result = simulator.run(qobj).result()
    counts = result.get_counts(qc)
    
    return counts

def refine_pi_estimate(previous_estimate, new_measurement, learning_rate=0.1):
    # Convert the measurement to phase
    measured_phase = int(new_measurement, 2) / (2 ** len(new_measurement))
    # Update the estimate of pi
    updated_pi = previous_estimate + learning_rate * (measured_phase * 4 - previous_estimate)
    return updated_pi

# Initialize variables
pi_estimate = 3.0  # Start with an arbitrary estimate of pi
estimates = [pi_estimate]
bit_precision = 5  # Number of bits for phase estimation

# Perform multiple iterations to refine the estimate
for _ in range(10):  # Number of iterations to refine the estimate
    counts = run_quantum_phase_estimation(bit_precision)
    # Take the most probable result
    most_probable_result = max(counts, key=counts.get)
    pi_estimate = refine_pi_estimate(pi_estimate, most_probable_result)
    estimates.append(pi_estimate)

# Plot the evolution of the pi estimate
plt.figure(figsize=(10, 5))
plt.plot(estimates, marker='o')
plt.title('Convergence of Pi Estimation Using Quantum Phase Estimation')
plt.xlabel('Iteration')
plt.ylabel('Estimate of Pi')
plt.grid(True)
plt.show()
