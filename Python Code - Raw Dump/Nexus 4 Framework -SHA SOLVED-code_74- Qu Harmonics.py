import hashlib

def mine_hash(data, initial_nonce, target_prefix="0000"):
    """
    Simulates a mining process using a nonce to find a hash that meets the target condition.

    Parameters:
        data (str): The input data (e.g., "HelloHello").
        initial_nonce (int): The starting value for the nonce.
        target_prefix (str): The target prefix (e.g., "0000" for Bitcoin-like difficulty).

    Returns:
        dict: Contains the resulting nonce, hash, and attempts.
    """
    nonce = initial_nonce
    max_nonce = 2 ** 32 - 1  # Max unsigned 32-bit integer
    attempts = 0

    while nonce <= max_nonce:
        # Combine data and nonce
        input_data = f"{data}{nonce}".encode('utf-8')
        
        # Generate the hash
        hash_result = hashlib.sha256(input_data).hexdigest()
        attempts += 1

        # Check if the hash meets the target condition
        if hash_result.startswith(target_prefix):
            return {
                "nonce": nonce,
                "hash": hash_result,
                "attempts": attempts
            }

        # Increment the nonce
        nonce += 1

    # If no valid nonce is found
    return {
        "nonce": None,
        "hash": None,
        "attempts": attempts
    }

# Test inputs
data = "HelloHello"
initial_nonce = -1481458232  # Starting nonce
result = mine_hash(data, initial_nonce, target_prefix="0000")

# Output results
print("Mining Results:")
print(f"Nonce: {result['nonce']}")
print(f"Hash: {result['hash']}")
print(f"Attempts: {result['attempts']}")
