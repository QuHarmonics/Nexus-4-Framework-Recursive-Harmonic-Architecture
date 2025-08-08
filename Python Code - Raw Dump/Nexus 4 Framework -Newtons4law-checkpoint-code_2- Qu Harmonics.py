import hashlib

def sha256_hex(text):
    return hashlib.sha256(text.encode()).hexdigest()

print("hellp:", sha256_hex("hellp"))
print("hello :", sha256_hex("hello"))