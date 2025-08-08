hex_to_bin = lambda x: bin(int(x, 16))[2:]
bin_to_int = lambda x: int(x, 2)

hash_hex = "8f434346648f6b96df89dda901c5176b10a6d83961dd3c1ac88b59b2dc327aa4"
hash_bin = hex_to_bin(hash_hex)
hash_int = bin_to_int(hash_bin)

print("Hash (hex):", hash_hex)
print("Hash (bin):", hash_bin)
print("Hash (int):", hash_int)