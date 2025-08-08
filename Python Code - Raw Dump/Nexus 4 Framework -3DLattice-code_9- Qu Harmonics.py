# Step 3: Enhanced Feedback Correction
def feedback_correction(lattice, binary_data, harmonic_constant=HARMONIC_CONSTANT):
    retrieved_data = retrieve_from_lattice(lattice, harmonic_constant, data_length=len(binary_data))
    error = (binary_data - retrieved_data) / 255.0  # Normalize the error
    
    for idx, value in enumerate(error):
        x, y, z = idx % lattice.shape[0], (idx // lattice.shape[0]) % lattice.shape[1], idx // (lattice.shape[0] ** 2)
        lattice[x, y, z] += value * harmonic_constant * 0.5  # Scale down the correction to avoid overshoot
    
    return lattice

# Step 4: Dynamically Adjust Reflective Gain
def apply_reflective_gain(lattice, gain_factor=0.1, iteration=1):
    center = lattice.shape[0] // 2
    dynamic_gain = gain_factor / (iteration + 1)  # Reduce gain with each iteration
    for x in range(lattice.shape[0]):
        for y in range(lattice.shape[1]):
            for z in range(lattice.shape[2]):
                distance = np.sqrt((x - center)**2 + (y - center)**2 + (z - center)**2)
                lattice[x, y, z] += dynamic_gain / (1 + distance)
    return lattice

# Step 5: Visualize Byte Differences
def visualize_byte_differences(original_data, retrieved_data):
    differences = original_data - retrieved_data
    plt.figure(figsize=(12, 6))
    plt.bar(range(len(differences)), differences, color='orange', alpha=0.7)
    plt.title("Byte-wise Differences Between Original and Retrieved Data")
    plt.xlabel("Byte Index")
    plt.ylabel("Difference")
    plt.show()

# Main Execution with Enhanced Logic
if __name__ == "__main__":
    # Load binary bits from your file
    with open(r'd:\\colecovision.rom', 'rb') as file:
        binary_data = np.frombuffer(file.read(), dtype=np.uint8)  # Read binary as bytes

    # Initialize the lattice
    lattice, data_length = initialize_lattice(binary_data)

    for iteration in range(1, MAX_ITERATIONS + 1):
        print(f"Iteration {iteration}")
        
        # Apply feedback correction
        lattice = feedback_correction(lattice, binary_data)
        
        # Apply dynamic reflective gain
        lattice = apply_reflective_gain(lattice, gain_factor=0.05, iteration=iteration)
        
        # Retrieve data
        retrieved_data = retrieve_from_lattice(lattice, data_length=data_length)

        # Visualize the lattice
        visualize_lattice(lattice, iteration)

        # Visualize byte-wise differences
        visualize_byte_differences(binary_data[:100], retrieved_data[:100])  # First 100 bytes for clarity

        # Outputs: Compare Original and Retrieved Data
        print("Lattice Shape:", lattice.shape)
        print("Original Data (First 10 Bytes):", binary_data[:10])
        print("Retrieved Data (First 10 Bytes):", retrieved_data[:10])

        # Validate if original and retrieved data match
        data_matches = np.array_equal(binary_data, retrieved_data)
        print("Data matches:", data_matches)
        if data_matches:
            print("Data successfully recovered!")
            break
        else:
            print("Differences detected in data.")
