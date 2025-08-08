# Hexadecimal to English Mapping
hex_to_text = {
    "41": "A", "42": "B", "43": "C", "44": "D",  # Standard ASCII for A-D
    "20": " ",                                   # Space
    "31": "1", "32": "2", "33": "3",             # Digits
}

def hex_to_english(hex_sequence):
    return ''.join([hex_to_text.get(h, '?') for h in hex_sequence.split()])

# Test
hex_sequence = "41 42 20 43 44"
print(hex_to_english(hex_sequence))  # Outputs: "AB CD"
