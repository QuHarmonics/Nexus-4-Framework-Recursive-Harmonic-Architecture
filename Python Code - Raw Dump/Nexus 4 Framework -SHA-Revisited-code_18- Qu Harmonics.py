import hashlib
import matplotlib.pyplot as plt
import numpy as np

# === SHA Pressure Tracker ===
def sha256_pressure(data):
    """
    Apply SHA-256 twice, return hash and pseudo-pressure from each byte delta.
    """
    first_hash = hashlib.sha256(data.encode()).digest()
    second_hash = hashlib.sha256(first_hash).digest()

    # Simulate "pressure" as deviation from average byte value
    avg = sum(second_hash) / len(second_hash)
    pressure = [abs(b - avg) for b in second_hash]
    return second_hash.hex(), pressure


# === Recursive Constraint Removal Simulation ===
def unfold_recursively(seed, rounds=20):
    results = []
    input_data = seed
    for i in range(rounds):
        hex_out, pressure = sha256_pressure(input_data)
        results.append({
            'round': i,
            'input': input_data,
            'hash': hex_out,
            'pressure': pressure
        })

        # Simulate release: fold hash back into text
        input_data = hex_out[:32]  # Take first 16 bytes (compressed loop)
    return results


# === Visualization ===
def plot_pressure(results):
    plt.figure(figsize=(10, 6))
    for result in results:
        plt.plot(result['pressure'], label=f"Round {result['round']}", alpha=0.6)
    plt.title("SHA-256 Byte Pressure per Round")
    plt.xlabel("Byte Index")
    plt.ylabel("Deviation from Mean (Pseudo-Pressure)")
    plt.legend()
    plt.grid(True)
    plt.show()


# === Run Simulation ===
if __name__ == "__main__":
    seed_phrase = "HELLO"
    results = unfold_recursively(seed_phrase, rounds=10)
    plot_pressure(results)

    for r in results[-3:]:  # Show last few rounds for reflection
        print(f"Round {r['round']}:\n Input: {r['input']}\n Hash: {r['hash'][:32]}...\n")
