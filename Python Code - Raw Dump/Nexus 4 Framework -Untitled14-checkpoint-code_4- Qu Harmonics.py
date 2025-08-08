def calculate_container_size(target, past_values):
    """
    Calculate the container size based on the target number and relationship with past values.
    
    Args:
    target (int): The target number.
    past_values (list): A list of past values.
    
    Returns:
    int: The container size.
    """
    # Calculate the container size based on the target number and past values
    container_size = max(target, max(past_values)) + len(past_values)
    return container_size


def populate_container(container_size, past_values):
    """
    Populate the container with valid combinations using recursive rules.
    
    Args:
    container_size (int): The container size.
    past_values (list): A list of past values.
    
    Returns:
    list: A list of valid combinations.
    """
    # Initialize an empty list to store valid combinations
    combinations = []
    
    # Define recursive rules to generate valid combinations
    def recursive_rules(current_value, remaining_values, depth):
        if depth == 0:
            # If the maximum depth is reached, add the current value to the combinations list
            if 0 <= current_value <= container_size:
                combinations.append(current_value)
        else:
            # Otherwise, recursively apply the rules to the remaining values
            for i in range(len(remaining_values)):
                for operator in [lambda x, y: x + y, lambda x, y: x - y, lambda x, y: x * y, lambda x, y: x / y]:
                    new_value = operator(current_value, remaining_values[i])
                    recursive_rules(new_value, remaining_values[:i] + remaining_values[i+1:], depth - 1)
    
    # Start the recursion with an initial value of 0, the past values, and a maximum depth
    recursive_rules(0, past_values, len(past_values))
    
    return combinations


# Example usage:
target = 10
past_values = [3, 1, 4, 1, 5, 9]

container_size = calculate_container_size(target, past_values)
combinations = populate_container(container_size, past_values)

print(f"Container size: {container_size}")
print(f"Combinations: {combinations}")