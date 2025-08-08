import hashlib
from mpmath import mp
import matplotlib.pyplot as plt

# Load 100,000 digits of π
mp.dps = 100_010
pi_digits = str(mp.pi)[2:100_010]

def extract_pi_byte(index: int) -> str:
    if index + 8 > len(pi_digits):
        return None
    return pi_digits[index:index+8]

def drift_and_echo(byte: str):
    deltas = [abs(int(byte[i+1]) - int(byte[i])) for i in range(7)]
    echo = ''.join([chr((d % 26) + 97) for d in deltas])
    avg_drift = sum(deltas) / len(deltas)
    sti = round(1 - avg_drift / 9, 3)
    return deltas, echo, sti

def run_recursive_kernel(seed: str, max_iter=25, sti_threshold=0.7):
    print(f"🔁 Running Recursive Kernel: {seed}")
    H = seed
    history = []

    for i in range(max_iter):
        nonce = f"N{i}"
        double_hash = hashlib.sha256(hashlib.sha256((H + nonce).encode()).digest()).hexdigest()
        index = int(double_hash[:6], 16) % (len(pi_digits) - 8)
        byte = extract_pi_byte(index)
        deltas, echo, sti = drift_and_echo(byte)

        print(f"🧠 Iter {i} | π@{index} → {byte} → {echo} | Δπ: {deltas} | STI: {sti}")
        history.append((i, index, byte, echo, deltas, sti))

        if sti >= sti_threshold:
            print(f"✅ ZPHC reached at iteration {i} with STI = {sti}")
            break

        H = double_hash

    return history

# Run it
history = run_recursive_kernel("RECURSE-ME-01")

# Optional: Plot STI curve
try:
    import matplotlib.pyplot as plt
    iters = [h[0] for h in history]
    stis = [h[5] for h in history]
    plt.plot(iters, stis, marker='o')
    plt.axhline(0.7, color='r', linestyle='--', label='ZPHC Threshold')
    plt.title("STI Over Recursive Iterations")
    plt.xlabel("Iteration")
    plt.ylabel("STI")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
except ImportError:
    pass
