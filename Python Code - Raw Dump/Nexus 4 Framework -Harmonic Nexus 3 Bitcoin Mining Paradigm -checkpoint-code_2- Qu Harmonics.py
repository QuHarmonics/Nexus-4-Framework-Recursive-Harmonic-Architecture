import hashlib
import matplotlib.pyplot as plt

# === Block Header Constants ===
version = bytes.fromhex('01000000')
prev_hash = bytes.fromhex('00' * 32)
merkle_root = bytes.fromhex('3BA3EDFD7A7B12B27AC72C3E67768F617FC81BC3888A51323A9FB8AA4B1E5E4A')[::-1]
time = bytes.fromhex('29ab5f49')
bits = bytes.fromhex('ffff001d')

# === Difficulty ===
TARGET_ZERO_BITS = 48  # Realistic modern mining target

# === Count leading zero bits ===
def count_leading_zero_bits(h_bytes):
    h_int = int.from_bytes(h_bytes, 'big')
    return 256 - h_int.bit_length() if h_int != 0 else 256

# === Harmonic Miner ===
def harmonic_mine(initial_nonce, max_iterations=1_000_000_000):
    nonce = initial_nonce
    best_zero_bits = 0
    history_best = []

    for iteration in range(max_iterations):
        nonce_bytes = nonce.to_bytes(4, 'little', signed=False)
        header = version + prev_hash + merkle_root + time + bits + nonce_bytes

        h1 = hashlib.sha256(header).digest()
        h2 = hashlib.sha256(h1).digest()

        zero_bits = count_leading_zero_bits(h2)
        if zero_bits > best_zero_bits:
            best_zero_bits = zero_bits
            print(f"[+] New Best: Nonce={nonce}, Zeros={zero_bits}, Hash={h2.hex()}")

        history_best.append(best_zero_bits)

        if zero_bits >= TARGET_ZERO_BITS:
            print(f"\n🎉 SUCCESS! Nonce={nonce} meets target with {zero_bits} zeros!")
            break

        # Harmonic Feedback Update
        feedback = int.from_bytes(h2[:4], 'little')
        nonce = nonce ^ feedback  # Phase-shift next input

    return history_best

# === Run Miner ===
initial_nonce = 3141592653 % (2**32)
history = harmonic_mine(initial_nonce)

# === Plot Progress ===
plt.plot(history)
plt.title("Harmonic Mining Progress")
plt.xlabel("Iterations")
plt.ylabel("Leading Zero Bits")
plt.grid(True)
plt.show()
