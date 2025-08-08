import math

def nexus2_hash(input_str):
    # 1. Convert input to its hexadecimal representation.
    pre_hash = input_str.encode('utf-8').hex()  
    # 2. Pad the hex string to a fixed length (for example, 64 hex digits)
    padded = pre_hash.ljust(64, '0')  # This is our "padded pre-hash"
    
    # 3. Split the padded string into two equal halves (A and B)
    A_hex = padded[:32]
    B_hex = padded[32:]
    # Convert these halves to integers
    A = int(A_hex, 16)
    B = int(B_hex, 16)
    
    # 4. Define constant parameters:
    C = 0.35               # Harmonic constant (universal attractor)
    LenC = len(str(C).replace('.', ''))  # e.g., "0.35" -> "35" gives 2 digits
    k = 0.1                # Feedback constant (for dynamic adjustments)
    
    # 5. Define a scaling factor AX (for demonstration, use the normalized average)
    # In practice, this would be derived from the kinetic energy of the system.
    AX = (A + B) / 2.0  
    
    # 6. Compute the exponential adjustment factor.
    exp_term = math.exp(-10 * (AX - C))
    
    # 7. Compute the Universal Formula value F.
    F_value = (A**2 + B**2) * LenC * (1 + exp_term)
    
    # 8. (Optional) Iterate the process to simulate recursive refinement:
    # For simplicity, we show a single-pass result.
    
    # Convert final numeric value to hexadecimal (or any normalized form)
    final_hash = hex(int(F_value))
    return final_hash

# Example usage:
print("Final hash for 'hello, how are you':", nexus2_hash("hello, how are you"))
