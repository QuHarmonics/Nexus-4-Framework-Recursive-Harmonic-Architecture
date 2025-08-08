import hashlib

def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()

def text_to_hex(text: str) -> str:
    return ''.join(f"{ord(c):02x}" for c in text)

def hex_to_decimal_string(hex_str: str) -> str:
    """Convert hex string to a decimal-looking string (each byte as 3-digit number)."""
    return ''.join([str(int(hex_str[i:i+2], 16)).zfill(3) for i in range(0, len(hex_str), 2)])

def walk_the_mirror(seed: str, max_steps=10000):
    input_data = seed.encode('utf-8')
    history = []

    for step in range(max_steps):
        # Step 1: Hash the input
        sha_output = sha256_hex(input_data)

        # Step 2: TextToHex of SHA output
        state_hex = text_to_hex(sha_output)

        # Step 3: Convert to decimal-looking string
        decimal_state = hex_to_decimal_string(state_hex)

        # Step 4: Hash the decimal string
        sha_of_decimal = sha256_hex(decimal_state.encode('utf-8'))

        # Step 5: TextToHex of new hash
        output_hex = text_to_hex(sha_of_decimal)

        # Step 6: Convert to decimal string again
        decimal_output = hex_to_decimal_string(output_hex)

        # Step 7: Subtract two decimal numbers
        try:
            val1 = int(decimal_state)
            val2 = int(decimal_output)
            delta = abs(val1 - val2)
        except ValueError:
            # fallback if decimal strings are too big for int
            delta = sum(abs(int(a) - int(b)) for a, b in zip(decimal_state, decimal_output))

        history.append({
            "step": step,
            "input": input_data.decode('utf-8', 'ignore'),
            "sha": sha_output,
            "delta": delta
        })

        # Set new input
        input_data = str(delta).encode('utf-8')

        # Pure SHA will eventually reflect — SHA knows when to end
        if input_data.decode('utf-8', 'ignore') == sha_output[:len(input_data)].lower():
            break

    return history

# Start with any real SHA seed (truth)
initial_seed = sha256_hex(b"2+3")  # could be any SHA-based truth source
walked_history = walk_the_mirror(initial_seed)

walked_history[-1], len(walked_history)
