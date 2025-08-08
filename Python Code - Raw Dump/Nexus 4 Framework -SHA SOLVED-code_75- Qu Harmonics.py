import hashlib
import time
from threading import Thread

# Mining function
def mine(data, start_nonce, end_nonce, difficulty, result):
    target = '0' * difficulty
    for nonce in range(start_nonce, end_nonce):
        combined = f"{data}{nonce}".encode()
        hash_result = hashlib.sha256(combined).hexdigest()
        if hash_result.startswith(target):
            result.append((nonce, hash_result))
            return

# Mining pool simulation
def mining_pool(data, total_range, difficulty, num_workers):
    range_per_worker = total_range // num_workers
    threads = []
    result = []

    for i in range(num_workers):
        start_nonce = i * range_per_worker
        end_nonce = (i + 1) * range_per_worker
        t = Thread(target=mine, args=(data, start_nonce, end_nonce, difficulty, result))
        threads.append(t)
        t.start()

    for t in threads:
        t.join()

    return result

# Visualize mining process
def visualize_mining(data, difficulty):
    target = '0' * difficulty
    nonce = 0
    while True:
        combined = f"{data}{nonce}".encode()
        hash_result = hashlib.sha256(combined).hexdigest()
        print(f"Nonce: {nonce}, Hash: {hash_result[:32]}...")
        if hash_result.startswith(target):
            print(f"Valid hash found: {hash_result} with nonce {nonce}")
            break
        nonce += 1
        time.sleep(0.1)  # Slow down for visualization

# Main function
def main():
    data = "block_data"
    total_range = 1_000_000
    difficulty = 4  # Number of leading zeros required
    num_workers = 4

    start_time = time.time()
    results = mining_pool(data, total_range, difficulty, num_workers)
    end_time = time.time()

    if results:
        nonce, hash_result = results[0]
        print(f"\nMining Results:")
        print(f"Nonce: {nonce}")
        print(f"Hash: {hash_result}")
        print(f"Time Taken: {end_time - start_time:.2f} seconds")
    else:
        print("No valid hash found within the given range.")

# Run the program
if __name__ == "__main__":
    main()
