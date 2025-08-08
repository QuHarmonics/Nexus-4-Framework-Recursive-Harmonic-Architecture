def simulate_assembly(input_data):
    # Initialize state
    registers = {
        "eax": 0,
        "ebx": 0,
        "ecx": 0,
        "edx": 0,
        "esi": 0,
        "edi": 0,
        "esp": 0,
        "ebp": 0,
        "eflags": 0,  # For flags like ZF, CF, etc.
    }
    memory = {}  # Simulated memory

    # Input setup
    input_bytes = [int(b) for b in input_data]  # e.g., [1, 8]
    memory[0] = input_bytes[0]  # Simulating some memory usage
    memory[1] = input_bytes[1]

    # Function to evaluate memory address
    def evaluate_address(address_str):
        """Evaluate memory address involving registers and offsets."""
        # Replace register names with their values
        for reg in registers:
            address_str = address_str.replace(reg, str(registers[reg]))
        # Evaluate the expression
        return eval(address_str)

    # Simulate instructions
    instructions = [
        {"op": "sbb", "dst": "[edx]", "src": "ah"},
        {"op": "sub", "dst": "eax", "src": "0x584d4238"},
        {"op": "bound", "dst": "ebp", "src": "[ebp+0x78]"},
        {"op": "or", "dst": "[ebp-0x47525d68]", "src": "0xc2"},
        {"op": "les", "dst": "edi", "src": "[ebp+esi*4-0x69615953]"},
        {"op": "pop", "dst": "[edi+0x69717880]"},
        {"op": "bound", "dst": "ebx", "src": "[edx+0x53]"},
        {"op": "dec", "dst": "esp"},
    ]

    # Simulate
    for idx, instr in enumerate(instructions):
        op = instr["op"]
        dst = instr["dst"]
        src = instr.get("src", None)  # Default to None if "src" key is missing

        try:
            if op == "sbb":
                reg = dst.strip("[]")  # Simplify for demonstration
                registers[reg] -= registers.get(src, 0)
                print(f"After {op} {dst},{src}: {reg}={registers[reg]}")

            elif op == "sub":
                if dst in registers:
                    registers[dst] -= int(src, 16)
                    print(f"After {op} {dst},{src}: {dst}={registers[dst]}")

            elif op == "dec":
                if dst in registers:
                    registers[dst] -= 1
                    print(f"After {op} {dst}: {dst}={registers[dst]}")

            elif op == "or":
                if "[" in dst:  # Memory
                    addr = evaluate_address(dst.strip("[]"))
                    memory[addr] = memory.get(addr, 0) | int(src, 16)
                    print(f"After {op} {dst},{src}: memory[{addr}]={memory[addr]}")

            elif op == "les":
                # Simplify: Assume loading pointer
                addr = evaluate_address(src.strip("[]"))
                registers[dst] = memory.get(addr, 0)
                print(f"After {op} {dst},{src}: {dst}={registers[dst]}")

            elif op == "pop":
                reg = dst.strip("[]")
                registers[reg] = registers["esp"]
                registers["esp"] += 4
                print(f"After {op} {dst}: {reg}={registers[reg]}")

        except Exception as e:
            print(f"Error during instruction {idx}: {op} {dst},{src} -> {e}")

    print("\nFinal State:")
    print("Registers:", registers)
    print("Memory:", memory)

# Test input
simulate_assembly(input_data=["1", "8"])
