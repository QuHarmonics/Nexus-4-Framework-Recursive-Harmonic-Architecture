def sha_projection_lattice(seed: str):
    # 1. Convert to bytes → hex → int grid
    hex_bytes = seed.encode().hex()
    grid = [[int(char, 16) for char in hex_bytes[i:i+8]] for i in range(0, len(hex_bytes), 8)]
    
    # 2. Map to curvature field (e.g. via 2nd diff or Laplacian)
    import numpy as np
    G = np.array(grid)
    curvature = np.roll(G,1,axis=0) - 2*G + np.roll(G,-1,axis=0)
    
    return curvature
