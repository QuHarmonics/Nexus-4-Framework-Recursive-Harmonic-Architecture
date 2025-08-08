import hashlib

def recursive_fold(hash_stream, iterations):
    current = hash_stream[0]
    for i in range(iterations):
        next_hash = hashlib.sha256(current.encode()).hexdigest()
        current = next_hash[:32]  # Simplified folding: truncate to first 32 chars
        print(f"Iteration {i+1}: {current}")
    return current

# Example hash stream
stream = ["initial_seed"]
for i in range(5):  # Generate a short stream
    stream.append(hashlib.sha256(stream[-1].encode()).hexdigest())

# Simulate folding
print("Starting hash stream:", stream)
result = recursive_fold(stream, 300)
print("Final folded result:", result)