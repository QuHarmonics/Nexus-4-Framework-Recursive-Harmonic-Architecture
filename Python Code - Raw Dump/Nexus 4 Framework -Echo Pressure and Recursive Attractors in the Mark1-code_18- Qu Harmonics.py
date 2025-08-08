def recursive_sha_square(seed, rounds=100):
    state = seed.encode()
    for i in range(rounds):
        # SHA
        h = hashlib.sha256(state).hexdigest()
        ascii_bytes = bytes([ord(c) for c in h])
        
        # Interpretations
        as_bin = ''.join(f'{b:08b}' for b in ascii_bytes)
        as_dec = [int(f'{b}', 16) for b in h]

        # Log the corners
        print(f"[{i}] SHA HEX: {h[:16]}... [END: {h[-2:]}]")
        print(f"   ▸ Bin head: {as_bin[:32]}...")
        print(f"   ▸ Dec tail: {as_dec[-4:]}")
        print(f"   ▸ Hex tail: {ascii_bytes[-4:].hex()}")

        # Feedback
        state = ascii_bytes  # You can rotate: binary, dec, etc.

    print("Terminus reached.")
