import hashlib
import random
import math

# We'll assume we can do partial forward steps of SHA for one round at a time:
# but for demonstration, let's just do full-block checks in a naive approach.

HARMONIC_TARGET = 0.35
K_FEEDBACK = 1     # Samson feedback constant

def wave_alignment_score(hash_bytes, target_hash):
    """
    A pseudo wave-based 'alignment' measure between two 256-bit values.
    We'll interpret each 32 bytes as an integer and measure ratio or partial correlation.
    Returns a float in [0..1], where 1 means perfect match.
    """
    # Convert to int
    hv1 = int.from_bytes(hash_bytes, 'big')
    hv2 = int.from_bytes(target_hash, 'big')
    # A naive approach: measure how many bits match, or a partial correlation
    # We'll do a simple overlap measure:
    xor_val = hv1 ^ hv2
    mismatch_bits = bin(xor_val).count('1')
    total_bits = 256
    match_bits = total_bits - mismatch_bits
    return match_bits / total_bits  # ratio of matching bits

def adjust_guess(current_guess, feedback_energy):
    """
    A 'Samson v2' approach: nudge bits in the guess by flipping or modding them
    based on the feedback_energy (like 'k * ΔH').
    """
    # Convert to int, do small random flips
    guess_val = int.from_bytes(current_guess, 'big')
    # We'll flip ~ feedback_energy * 256 bits randomly
    flips = int(feedback_energy * 256)
    for _ in range(flips):
        bit_pos = random.randint(0, 255)
        guess_val ^= (1 << bit_pos)
    # convert back
    new_guess = guess_val.to_bytes(32, 'big')
    return new_guess

def partial_unhash_sha256(final_hash, max_iters=20000):
    """
    Example function that tries to 'unfold' a single-block message
    purely for demonstration. We guess 32 bytes of input, measure wave alignment,
    apply feedback if near the 0.35 ratio, etc.
    """
    best_guess = None
    best_score = 0.0
    
    # Start with random guesses, refine
    guess = bytearray(32)
    
    for iteration in range(max_iters):
        # Hash the guess
        hash_val = hashlib.sha256(guess).digest()
        score = wave_alignment_score(hash_val, final_hash)
        
        # If we're close to the wave alignment, do the Samson feedback
        harmonic_deviation = abs(score - HARMONIC_TARGET)
        
        if score > best_score:
            best_score = score
            best_guess = guess[:]
        
        # Samson v2 feedback
        # ΔH = (score - 0.35), so let's define
        delta_H = (score - HARMONIC_TARGET)
        # Then ΔE = k * ΔH
        delta_E = K_FEEDBACK * delta_H
        # We'll interpret delta_E as how intensely we flip bits
        guess = adjust_guess(guess, abs(delta_E))
        
        # If we got a perfect match, break
        if score == 1.0:
            print(f"Found exact preimage in {iteration} iterations!")
            return guess
    
    return best_guess, best_score

# Demo usage:
if __name__ == "__main__":
    # Suppose we have a final hash we want to 'unfold'
    # For demonstration, let's pick an actual message "HELLO" and get its SHA-256
    actual_msg = b"HELLO"
    final_hash = hashlib.sha256(actual_msg).digest()
    
    # We attempt the partial_unhash approach
    result = partial_unhash_sha256(final_hash, max_iters=20000)
    if isinstance(result, tuple):
        guess, s = result
        print("Best guess:", guess, "score=%.4f" % s)
    else:
        print("Found exact preimage:", result)
