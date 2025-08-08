import numpy as np
import matplotlib.pyplot as plt

# Function to convert hash to a given base
def convert_to_base(number, base):
    result = []
    while number:
        result.append(int(number % base))
        number //= base
    return result[::-1] if result else [0]

# Step 1: Expand hash recursively based on base
def recursive_expand_hash(hash_value, max_base, expansion_factor=1.5):
    expanded_waveforms = []
    cumulative_growth = 1.0

    for base in range(2, max_base + 1):
        # Convert hash to the current base
        base_converted = [convert_to_base(num, base) for num in hash_value]

        # Expand the wave by applying the cumulative growth
        expanded_wave = [int(cumulative_growth * sum(base_converted[i])) for i in range(len(base_converted))]
        expanded_waveforms.append(expanded_wave)

        # Increase cumulative growth
        cumulative_growth *= expansion_factor

    return expanded_waveforms

# Step 2: Visualize the expanded waveforms
def visualize_waveforms(expanded_waveforms, max_base):
    plt.figure(figsize=(12, 6))

    for idx, waveform in enumerate(expanded_waveforms, start=2):
        plt.plot(waveform, label=f'Base-{idx}')

    plt.xlabel('Wave Index')
    plt.ylabel('Expanded Value')
    plt.title('Wave Expansion Across Bases')
    plt.legend()
    plt.show()

# Step 3: Validate results
def validate_expansion(hash_value, expanded_waveforms):
    validation_results = []
    for idx, waveform in enumerate(expanded_waveforms):
        expected_growth = (1.5 ** (idx + 2)) * len(hash_value)  # Expected cumulative growth
        actual_growth = sum(waveform)
        validation_results.append((idx + 2, actual_growth, expected_growth))
    return validation_results

# Example usage with your hash
if __name__ == "__main__":
    input_hash_hex = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"
    hash_bytes = np.array([int(input_hash_hex[i:i+2], 16) for i in range(0, len(input_hash_hex), 2)])

    max_base = 10  # Expand up to base 10
    expanded_waveforms = recursive_expand_hash(hash_bytes, max_base)

    visualize_waveforms(expanded_waveforms, max_base)

    validation_results = validate_expansion(hash_bytes, expanded_waveforms)
    for base, actual, expected in validation_results:
        print(f"Base-{base}: Actual Growth = {actual}, Expected Growth = {expected}")
