import hashlib
import time
from threading import Thread
from queue import Queue

# Mining function
def mine(data, start_nonce, end_nonce, difficulty, result_queue):
    target = '0' * difficulty
    for nonce in range(start_nonce, end_nonce):
        combined = f"{data}{nonce}".encode()
        hash_result = hashlib.sha256(combined).hexdigest()
        if hash_result.startswith(target):
            result_queue.put((nonce, hash_result))
            return

# Mining pool simulation
def mining_pool(data, total_range, difficulty, num_workers):
    range_per_worker = total_range // num_workers
    threads = []
    result_queue = Queue()

    for i in range(num_workers):
        start_nonce = i * range_per_worker
        end_nonce = (i + 1) * range_per_worker
        t = Thread(target=mine, args=(data, start_nonce, end_nonce, difficulty, result_queue))
        threads.append(t)
        t.start()

    for t in threads:
        t.join()

    return list(result_queue.queue)

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

# Analyze hash assembly-level representation
def analyze_asm(hash_data):
    print("\nAnalyzing Hash as Assembly-Level Instructions...")
    asm_instructions = []
    try:
        for i in range(0, len(hash_data), 2):
            byte = hash_data[i:i+2]
            asm_instructions.append(f"0x{byte}")
    except Exception as e:
        print(f"Error processing assembly: {e}")
    print("Interpreted Assembly Instructions:")
    print(" ".join(asm_instructions))
    return asm_instructions

# Main function
def main():
    data = "block_data"
    total_range = 1_000_000
    difficulties = [4, 5, 6]  # Dynamic difficulties to test
    num_workers = 4

    for difficulty in difficulties:
        print(f"\n--- Mining with Difficulty {difficulty} ---")
        start_time = time.time()
        results = mining_pool(data, total_range, difficulty, num_workers)
        end_time = time.time()

        if results:
            nonce, hash_result = results[0]
            print(f"\nMining Results:")
            print(f"Nonce: {nonce}")
            print(f"Hash: {hash_result}")
            print(f"Time Taken: {end_time - start_time:.2f} seconds")
            print("Visualizing Mining Process...")
            visualize_mining(data, difficulty)
            print("\nAnalyzing Hash Assembly Representation...")
            analyze_asm(hash_result)
        else:
            print("No valid hash found within the given range.")

# Run the program
if __name__ == "__main__":
    main()
