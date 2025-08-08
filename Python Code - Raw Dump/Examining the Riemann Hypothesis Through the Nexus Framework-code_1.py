import hashlib

def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def hex_to_text(hex_str: str) -> str:
    """Convert hex string to text."""
    bytes_data = bytes.fromhex(hex_str)
    try:
        return bytes_data.decode('utf-8', errors='ignore')
    except UnicodeDecodeError:
        return ""

def hex_to_decimal_string(hex_str: str) -> str:
    """Convert hex string to a decimal-only representation (numeric string)."""
    return ''.join([str(int(hex_str[i:i+2], 16)).zfill(3) for i in range(0, len(hex_str), 2)])

def extract_numeric_values(decimal_string: str):
    """Extract meaningful numeric segments from a long decimal string."""
    if len(decimal_string) < 6:
        return 0, 0
    mid = len(decimal_string) // 2
    return int(decimal_string[:mid]), int(decimal_string[mid:])

def recursive_sha_feedback(seed: str, max_steps=100):
    input_data = seed.encode('utf-8')
    history = []

    for _ in range(max_steps):
        # Step 1: Hash the input
        sha_output = sha256_hex(input_data)

        # Step 2: Treat hash as hex → decode as text
        interpreted_text = hex_to_text(sha_output)

        # Step 3: Convert the original hex output as decimal-only string
        numeric_str = hex_to_decimal_string(sha_output)

        # Step 4: Extract two parts, subtract, send difference back as input
        val1, val2 = extract_numeric_values(numeric_str)
        delta = abs(val1 - val2)

        history.append((input_data.decode('utf-8', 'ignore'), sha_output, delta))

        # Prepare next input
        input_data = str(delta).encode('utf-8')

        # Optional: Check if input equals output in any way
        if input_data.decode('utf-8', 'ignore') == sha_output[:len(input_data)].lower():
            break

    return history

# Try it with a seed string
seed_value = "2+3"
result_history = recursive_sha_feedback(seed_value)

result_history[-1], len(result_history)
