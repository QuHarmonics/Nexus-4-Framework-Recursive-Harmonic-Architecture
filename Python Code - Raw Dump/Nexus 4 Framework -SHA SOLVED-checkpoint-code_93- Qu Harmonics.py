# Re-attempt base conversion after resolving execution state reset issue

def convert_to_bases(number):
    return {
        2: bin(number)[2:], 3: format(number, 'o'), 4: format(number, 'x'),
        5: format(number, 'b'), 6: format(number, 'o'), 7: format(number, 'x'),
        8: oct(number)[2:], 9: format(number, 'o'), 10: str(number),
        11: format(number, 'x'), 12: format(number, 'o'), 16: hex(number)[2:]
    }

# Convert numbers 1 and 4 to different bases
num_1_bases = convert_to_bases(1)
num_4_bases = convert_to_bases(4)

num_1_bases, num_4_bases
