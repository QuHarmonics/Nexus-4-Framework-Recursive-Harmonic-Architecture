# Function to identify and fill gaps, highlighting added nucleotides
def identify_and_fill_gaps_with_highlight(sequence, expected_length):
    """
    Identifies gaps in the sequence based on missing nucleotides, fills them symmetrically,
    and highlights the added nucleotides.
    
    Args:
    - sequence (str): The original nucleotide sequence.
    - expected_length (int): The expected length of the sequence after filling gaps.
    
    Returns:
    - filled_sequence (str): The sequence after filling the missing nucleotides.
    - highlighted_sequence (str): The sequence with added nucleotides highlighted.
    - total_pairs (int): Total number of pairs in the filled sequence.
    - added_pairs (int): Number of pairs formed by the inserted nucleotides.
    """
    original_length = len(sequence)
    missing_length = expected_length - original_length

    if missing_length <= 0:
        return sequence, sequence, original_length // 2, 0

    # Analyze sequence for gaps
    gaps = []
    for i in range(1, len(sequence)):
        if ord(sequence[i]) - ord(sequence[i - 1]) > 1:  # Gap detected
            gaps.append((i, ord(sequence[i]) - ord(sequence[i - 1]) - 1))

    # Redistribute missing nucleotides symmetrically into the sequence
    filled_sequence = list(sequence)
    highlighted_sequence = list(sequence)  # For storing highlighted sequence
    remaining_to_fill = missing_length

    for gap_index, gap_size in gaps:
        fill_size = min(gap_size, remaining_to_fill)
        filler = 'N' * fill_size  # Use 'N' for missing nucleotides
        remaining_to_fill -= fill_size

        # Insert missing nucleotides symmetrically
        filled_sequence.insert(gap_index, filler)
        highlighted_sequence.insert(gap_index, f"[{filler}]")  # Highlight added bases

        if remaining_to_fill <= 0:
            break

    # If there's still missing data after filling gaps, pad the end
    if remaining_to_fill > 0:
        filler = 'N' * remaining_to_fill
        filled_sequence.append(filler)
        highlighted_sequence.append(f"[{filler}]")  # Highlight added padding

    # Combine into complete sequences
    filled_sequence = ''.join(filled_sequence)
    highlighted_sequence = ''.join(highlighted_sequence)

    # Calculate pairs
    total_pairs = len(filled_sequence) // 2
    added_pairs = missing_length // 2

    return filled_sequence, highlighted_sequence, total_pairs, added_pairs
    
# Full Input Sequence
original_sequence = (
    "ATGGGGACGGAAGACTGCGATCACGAAGGGCGGTCGGTTGCGGCTCCCGTGGAGGTTACGGCGCTGTATG"
    "CGACCGACGGGTGCGTTATCACCTCCTCGCTCGCCCTCCTCACAAACTGCCTGCTGGGGGCCGAGCCGTT"
    "GTATATATTCAGCTACGACGCGTACCGGCCCGATGCGCCCAATGGCCCCACGGGCGCGCCCACCGAACAG"
    "GAGAGGTTCGAGGGGAGCCGGGCGCTCTACCGGGATGCGGGGGGGCTAAATGGCGATTCATTTCGGGTGA"
    "CCTTTTGTTTATTGGGGACGGAAGTGGGCGTGACCCACCACCCGAAAGGGCGCACCCGGCCCATGTTTGT"
    "GTGCCGCTTCGAGCGAGCGGACGACGTCGCCGTGCTCCAAGACGCCCTGGGCCGCGGGACCCCATTGCTC"
    "CCGGCCCACATCACAGCAACTCTGGACTTGGAGGCGACGTTTGCGCTCCACGCTAACATCATCATGGCTC"
    "TCACCGTGGCCATCGTCCACAACGCCCCCGCCCGCATCGGCAGCGGCAGCACCGCCCCCCTGTATGAGCC"
    "CGGCGAATCGATGCGCTCGGTCGTCGGGCGCATGTCCCTGGGGCAGCGCGGCCTCACCACGCTGTTCGTG"
    "CACCACGAGGCGCGCGTGCTGGCGGCGTACCGCCGGGCGTATTATGGGAGCGCCCAAAGCCCCTTTTGGT"
    "TTCTGAGCAAATTCGGCCCGGACGAAAAGAGCCTGGTGCTGGCCGCTAGGTACTACCTACTCCAGGCTCC"
    "GCGCTTGGGGGGCGCCGGAGCCACGTACGATCTGCAGGCCGTGAAAGACATCTGCGCGACCTACGCGATC"
    "CCCCACGACCCACGCCCCGACACCCTCAGTGCCGCGTCCTTGACCTCGTTCGCCGCCATCACTCGGTTCT"
    "GTTGCACGAGCCAGTACTCCCGCGGGGCCGCGGCCGCTGGGTTTCCGCTGTATGTGGAGCGCCGCATCGC"
    "CGCCGACGTACGCGAGACCGGCGCGCTGGAGAAGTTCATCGCCCACGATCGCAGCTGCCTGCGCGTGTCC"
    "GACCGGGAATTCATTACGTACATCTACCTGGCCCACTTTGAGTGCTTCAGCCCCCCGCGCCTGGCCACGC"
    "ATCTCCGGGCCGTGACCACCCACGACCCCAGCCCCGCGGCCAGCACGGAGCAGCCCTCGCCCCTGGGTCG"
    "GGAGGCGGTGGAACAGTTCTTCCGGCACGTGCGCGCCCAGCTGAACATCCGCGAGTACGTAAAGCAAAAC"
    "GTCACCCCCAGGGAAACCGCCCTGGCGGGAGACGCGGCCGCCGCCTACCTGCGCGCGCGCACGTATGCCC"
    "CGGCGGCCCTCACGCCCGCCCCCGCGTACTGCGGGGTCGCAGACTCGTCCACCAAAATGATGGGACGTCT"
    "GGCGGAAGCAGAAAGGCTCCTAGTCCCCCACGGCTGGCCCGCGTTCGCACCAACAACCCCCGGGGACGAC"
    "GCGGGGGGCGGCACTGCCGCCCCCCAGACCTGCGGAATCGTCAAGCGCCTCCTCAAGCTGGCCGCCACGG"
    "AGCAGCAGGGCACGACGCCCCCGGCGATCGCGGCTCTCATGCAGGACGCGTCGGTCCAAACCCCCCTGCC"
    "CGTGTACAGGATTACCATGTCCCCGACCGGCCAGGCGTTTGCCGCGGCGGCGCGGGACGACTGGGCCCGC"
    "GTGACGCGGGACGCGCGCCCGCCGGAAGCGACCGTGGTCGCGGACGCGGCGGCGGCGCCCGAGCCCGGCG"
    "CGCTCGGCCGGCGGCTCACGCGCCGCATTTGCGCCCGGGGCCCCGCGCTCCCCCCGGGCGGCCTGGCCGT"
    "CGGGGGCCAGATGTACGTGAACCGCAACGAGATCTTCAACGCCGCGCTGGCCGTTACGAACATCATCCTG"
    "GATCTGGACATCGCCCTGAAGGAGCCCGTCCCCTTTCCCCGGCTCCACGAGGCCCTGGGTCACTTTAGGC"
    "GCGGGGCGCTGGCGGCGGTTCAGCTGTTGTTTCCCGCGGCCCGCGTAGACCCCGACGCCTATCCCTGTTA"
    "TTTTTTCAAAAGCGCCTGTCGGCCCCGCGCGCCGCCCGTCTGTGCGGGCGACGGGCCCCTGGCCGGTGGC"
    "GACGACGGCGACGGGGACTGGTTCCCCGACGCCGGTGGTCCCGGCGACGAGGAGTGGGAGGAGGACACGG"
    "ACCCCATGGACACGACCCACGGCCCCCTCCCGGACGACGAGGCCGCGTACCTCGACCTGCTACACGAACA"
    "GATACCAGCGGCGACGCCCAGCGAACCGGACTCCGTCGTGTGTTCCTGCGCCGACAAGATCGGGCTGCGC"
    "GTGTGCCTACCGGTCCCCGCCCCGTACGTTGTGCACGGCTCCCTGACGATGCGTGGGGTGGCGAGGGTGA"
    "TCCAGCAGGCGGTGCTGTTGGACCGCGACTTCGTGGAGGCCGTAGGGAGCCACGTAAAGAACTTTTTGCT"
    "GATCGATACGGGCGTGTACGCCCACGGCCACAGCCTGCGCTTGCCGTATTTCGCCAAGATCGGCCCCGAC"
    "GGCTCCGCGTGCGGCCGGTTATTGCCCGTCTTCGTGATCCCCCCCGCGTGCGAGGACGTTCCGGCGTTCG"
    "TCGCCGCGCACGCCGACCCGCGGCGCTTCCACTTTCACGCCCCGCCCATGTTTTCCGCGGCCCCGCGGGA"
    "GATCCGCGTCCTCCACAGCCTGGGCGGGGACTATGTCAGCTTTTTCGAGAAGAAGGCGTCGCGCAACGCC"
    "CTGGAGCACTTTGGGCGACGCGAGACCCTGACGGAGGTTCTGGGCCGCTACGATGTGCGGCCCGACGCCG"
    "GGGAGACCGTGGAGGGGTTCGCGTCAGAACTGCTGGGGCGAATAGTCGCGTGCATCGAGGCCCACTTTCC"
    "CGAGCACGCGCGGGAATATCAGGCCGTGTCCGTTCGCCGGGCCGTCATTAAGGACGACTGGGTCCTGCTG"
    "CAGCTGATCCCCGGCCGCGGCGCCCTGAACCAAAGCCTCTCGTGTCTGCGCTTCAAGCACGGCAGGGCAA"
    "GTCGCGCGACGGCCCGGACCTTTCTCGCGCTGAGCGTCGGGACCAACAACCGCCTATGCGCGTCCCTGTG"
    "TCAGCAGTGCTTTGCCACTAAATGCGATAACAACCGCCTGCACACGCTGTTTACCGTCGATGCGGGCACG"
    "CCATGCTCGCGGTCCGCTCCCTCCAGCACCTCACGACCGTCATCTTCATAA")
    
expected_sequence_length = 4096

# Run the function
filled_sequence, highlighted_sequence, total_pairs, added_pairs = identify_and_fill_gaps_with_highlight(
    original_sequence, expected_sequence_length
)

# Display results
print("Original Sequence Length:", len(original_sequence))
print("Expected Sequence Length:", expected_sequence_length)
print("Filled Sequence Length:", len(filled_sequence))
print("Total Pairs in Filled Sequence:", total_pairs)
print("Added Pairs:", added_pairs)
print("\nHighlighted Sequence (First 500 bases):\n", highlighted_sequence)

# Optionally save the highlighted sequence to a file for further review
with open("highlighted_sequence.txt", "w") as file:
    file.write(highlighted_sequence)
