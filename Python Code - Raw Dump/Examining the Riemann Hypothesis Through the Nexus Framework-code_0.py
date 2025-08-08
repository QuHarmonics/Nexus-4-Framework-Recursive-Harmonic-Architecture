import hashlib
import random
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

# Initialize a 64-byte (512-bit) input seed
input_bytes = bytearray(random.getrandbits(8) for _ in range(64))

def sha256_hash(data: bytes) -> bytes:
    return hashlib.sha256(data).digest()

def score_alignment(input_bytes, hash_bytes, locked_indices):
    """Calculate how many mirrored bits match between input and SHA output, excluding locked."""
    match_score = 0
    for i in range(min(len(hash_bytes), len(input_bytes))):
        if i in locked_indices:
            continue
        mirrored_index = len(hash_bytes) - 1 - i
        if input_bytes[i] == hash_bytes[mirrored_index]:
            match_score += 1
    return match_score

def harmonize(input_bytes, max_iterations=500):
    locked = set()
    scores = []

    for step in range(max_iterations):
        original_input = deepcopy(input_bytes)
        hash_bytes = sha256_hash(bytes(input_bytes))

        # Score center to edge
        center = len(hash_bytes) // 2
        indices = list(range(center)) + list(range(center, len(hash_bytes)))

        for i in indices:
            if i in locked:
                continue

            mirror_index = len(hash_bytes) - 1 - i
            if input_bytes[i] != hash_bytes[mirror_index]:
                # Adjust input to match mirrored hash bit
                input_bytes[i] = hash_bytes[mirror_index]
            else:
                locked.add(i)  # Lock once matched

        current_score = score_alignment(input_bytes, hash_bytes, set())
        scores.append(current_score)

        # Break early if perfect match
        if current_score == len(hash_bytes):
            break

    return input_bytes, sha256_hash(bytes(input_bytes)), scores

# Run the harmonizer
final_input, final_hash, score_progression = harmonize(input_bytes)

# Plot the score progression
plt.figure(figsize=(10, 4))
plt.plot(score_progression, label='Alignment Score')
plt.xlabel('Iteration')
plt.ylabel('Matched Bits')
plt.title('SHA Reflection Harmony Progression')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

final_input.hex(), final_hash.hex()
