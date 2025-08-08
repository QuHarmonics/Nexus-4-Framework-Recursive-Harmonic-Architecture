print("Raw Text:", sha_input)
print("Hex-Decoded:", bytes.fromhex(hex_input))
print("Are they equal?", sha_input == bytes.fromhex(hex_input))