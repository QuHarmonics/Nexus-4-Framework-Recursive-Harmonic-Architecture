import hashlib
import numpy as np
import matplotlib.pyplot as plt

# Text encoding function
def text_to_prism_view(text):
    """Simulate refracting text through a prism by mapping characters to RGB."""
    prism_data = []
    for char in text:
        ascii_value = ord(char)
        r = (ascii_value * 3) % 256
        g = (ascii_value * 7) % 256
        b = (ascii_value * 11) % 256
        prism_data.append((r, g, b))
    return np.array(prism_data)

# Hash transformation function
def hash_text(text):
    """Hash text using SHA-256."""
    return hashlib.sha256(text.encode()).hexdigest()

# Reverse prism view function
def reverse_prism(prism_data):
    """Reverse the prism transformation to approximate the original text."""
    reconstructed_text = []
    for r, g, b in prism_data:
        # Approximation: reverse mapping RGB to ASCII values
        ascii_value = (r // 3 + g // 7 + b // 11) // 3
        reconstructed_text.append(chr(ascii_value % 128))  # Keep within printable ASCII
    return ''.join(reconstructed_text)

# Visualization function
def visualize_prism(prism_data, title):
    """Visualize the prism data as an image."""
    plt.imshow([prism_data], aspect='auto')
    plt.title(title)
    plt.axis('off')
    plt.show()

# Main simulation
def simulate_prism_hash():
    input_text = "Hello, world!"
    print("Original Text:", input_text)

    # Step 1: Encode text through the prism
    prism_data = text_to_prism_view(input_text)
    visualize_prism(prism_data, "Prism View of Text")

    # Step 2: Hash the original text
    hash_result = hash_text(input_text)
    print("Hash Result:", hash_result)

    # Step 3: Reverse the prism view
    reconstructed_text = reverse_prism(prism_data)
    print("Reconstructed Text (Approximation):", reconstructed_text)

# Run the simulation
simulate_prism_hash()
