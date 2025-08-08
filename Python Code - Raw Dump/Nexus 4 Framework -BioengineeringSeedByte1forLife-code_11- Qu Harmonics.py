from capstone import *

# Initialize Capstone for x86 in 32-bit mode
md = Cs(CS_ARCH_X86, CS_MODE_32)

def main():
    print("NASM Disassembler")
    print("Enter raw hex bytes (e.g., '9f803d3dde5f5f6a7b'):")
    
    # Get user input for raw hex bytes
    raw_hex = input("> ").strip()
    
    try:
        # Convert input into a bytes object
        code = bytes.fromhex(raw_hex)
    except ValueError:
        print("Invalid hex string. Please ensure it's valid hexadecimal characters.")
        return

    print("\nDisassembly Output:\n")
    
    # Disassemble the raw opcodes
    for i in md.disasm(code, 0x1000):  # 0x1000 is the default starting address
        print(f"{i.address:08x} {i.mnemonic} {i.op_str}")

if __name__ == "__main__":
    main()
