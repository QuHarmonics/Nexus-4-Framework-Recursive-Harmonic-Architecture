import math

def fold_with_zeta(binary, block_size):
    """
    Fold binary data using Zeta-axis alignment and triangular relationships.
    """
    zeta = []  # Zeta axis to align data
    folded_data = []
    metadata = []  # Store fold angles (e.g., theta) for reconstruction
    
    while len(binary) >= block_size:
        block = binary[:block_size]
        binary = binary[block_size:]  # Remove processed block
        
        # Calculate "triangle" relationships
        a = len(zeta)
        b = len(block)
        c = math.sqrt(a**2 + b**2)
        theta = math.degrees(math.atan(b / a)) if a > 0 else 90  # Angle in degrees
        
        # Fold the block into Zeta
        zeta.extend(block)
        folded_data.append((c, theta))  # Store diagonal and angle as metadata
    
    # Add remaining data
    if binary:
        zeta.extend(binary)
        folded_data.append((len(binary), 0))  # No angle for leftover data
    
    return folded_data, zeta

def unfold_with_zeta(folded_data, block_size):
    """
    Unfold binary data using Zeta-axis alignment and stored triangular relationships.
    """
    unfolded_data = []
    for c, theta in folded_data:
        if theta > 0:
            # Reconstruct the block using c and theta
            b = int(c * math.sin(math.radians(theta)))
            block = [1] * b  # Simulate block data (this needs refinement for actual data)
            unfolded_data.extend(block)
        else:
            # Remaining data
            unfolded_data.extend([1] * int(c))
    
    return unfolded_data

# Example binary string
binary_data = [1, 0, 1, 0, 1, 1, 0, 1] * 4
block_size = 8

# Fold data
folded, zeta = fold_with_zeta(binary_data, block_size)
print("Folded Data:", folded)
print("Zeta Axis:", zeta)

# Unfold data
unfolded = unfold_with_zeta(folded, block_size)
print("Unfolded Matches Original:", unfolded == binary_data)
