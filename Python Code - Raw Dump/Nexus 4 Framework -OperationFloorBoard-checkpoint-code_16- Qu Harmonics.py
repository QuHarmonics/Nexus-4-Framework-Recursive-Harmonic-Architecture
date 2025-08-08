# Import libraries for plotting and math
import math
import matplotlib.pyplot as plt

# Function to process input text
def process_input(input_text, base=16):
    x_values = []  # Gain (ratio of change for each bit length)
    y_values = []  # Frequency (difference between initial and iterative bit lengths)
    z_values = []  # Pythagorean result (a^2 + b^2 = c^2)

    initial_bit_length = len(bin(ord(input_text[0]))[2:])  # Bit length of the first character
    for i, char in enumerate(input_text):
        ascii_value = ord(char)  # ASCII value
        quotient = ascii_value // base  # Quotient in hex conversion
        remainder = ascii_value % base  # Remainder in hex conversion
        hex_value = f"{quotient:X}{remainder:X}"  # Hexadecimal representation

        # Gain: Ratio of change for bit lengths (current bit length vs previous bit length)
        current_bit_length = len(bin(ascii_value)[2:])  # Current bit length
        if i > 0:
            previous_bit_length = len(bin(ord(input_text[i - 1]))[2:])
            gain = current_bit_length / previous_bit_length if previous_bit_length != 0 else 0
        else:
            gain = 1  # No gain for the first character

        # Frequency: Difference between initial bit length and current bit length
        frequency = abs(current_bit_length - initial_bit_length)

        # Pythagorean result: sqrt(a^2 + b^2)
        a = quotient
        b = remainder
        c = math.sqrt(a**2 + b**2)

        # Append results
        x_values.append(gain)
        y_values.append(frequency)
        z_values.append(c)
    
    return x_values, y_values, z_values

# Input texts
inputs = ["2+8=11", "3+12=15", "1+1=6", "5*6=30"]  # Replace or extend as needed

# Create plots
plt.figure(figsize=(18, 12))

for idx, input_text in enumerate(inputs):
    x_values, y_values, z_values = process_input(input_text)

    # Plot Gain
    plt.subplot(4, 3, idx * 3 + 1)
    plt.plot(x_values, marker='o', label="Gain")
    plt.title(f"Gain for '{input_text}' (Bit Length Ratio Change)")
    plt.xlabel("Character Index")
    plt.ylabel("Gain")
    plt.grid()
    plt.legend()

    # Plot Frequency
    plt.subplot(4, 3, idx * 3 + 2)
    plt.plot(y_values, marker='o', label="Frequency", color='orange')
    plt.title(f"Frequency for '{input_text}' (Difference in Bit Lengths)")
    plt.xlabel("Character Index")
    plt.ylabel("Frequency")
    plt.grid()
    plt.legend()

    # Plot Transition Mix
    plt.subplot(4, 3, idx * 3 + 3)
    plt.plot(z_values, marker='o', label="Transition Mix", color='green')
    plt.title(f"Transition Mix for '{input_text}' (a^2 + b^2 = c^2)")
    plt.xlabel("Character Index")
    plt.ylabel("Transition Mix")
    plt.grid()
    plt.legend()

plt.tight_layout()
plt.show()
