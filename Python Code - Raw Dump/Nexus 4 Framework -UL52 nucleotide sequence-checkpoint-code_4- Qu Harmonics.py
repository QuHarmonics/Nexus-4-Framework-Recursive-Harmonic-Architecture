# Refined function to address gaps properly and minimize padding while grouping the sequence into pairs.

def identify_and_fill_gaps_properly(sequence, expected_length):
    """
    Refined version to detect gaps, fill them symmetrically, and minimize padding, 
    with the output grouped in pairs for readability.

    Args:
    - sequence (str): The original nucleotide sequence.
    - expected_length (int): The expected length of the sequence after filling gaps.

    Returns:
    - grouped_sequence (str): The sequence grouped in pairs for readability.
    - filled_sequence (str): The sequence after filling the missing nucleotides.
    - total_pairs (int): Total number of pairs in the filled sequence.
    - added_pairs (int): Number of pairs formed by the inserted nucleotides.
    """
    original_length = len(sequence)
    missing_length = expected_length - original_length

    if missing_length <= 0:
        grouped_sequence = " ".join([sequence[i:i+2] for i in range(0, len(sequence), 2)])
        return grouped_sequence, sequence, original_length // 2, 0

    # Detect gaps and distribute missing nucleotides
    gaps = []
    for i in range(1, len(sequence)):
        if ord(sequence[i]) - ord(sequence[i - 1]) > 1:  # Detect gaps
            gaps.append((i, ord(sequence[i]) - ord(sequence[i - 1]) - 1))

    filled_sequence = list(sequence)
    remaining_to_fill = missing_length

    # Fill gaps symmetrically
    for gap_index, gap_size in gaps:
        if remaining_to_fill <= 0:
            break
        fill_size = min(gap_size, remaining_to_fill)
        filler = 'N' * fill_size  # Use 'N' for missing nucleotides
        remaining_to_fill -= fill_size
        filled_sequence.insert(gap_index, filler)

    # If any nucleotides are still missing, add them to the end
    if remaining_to_fill > 0:
        filled_sequence.append('N' * remaining_to_fill)

    # Combine into complete sequence
    filled_sequence = ''.join(filled_sequence)

    # Group into pairs
    grouped_sequence = " ".join([filled_sequence[i:i+2] for i in range(0, len(filled_sequence), 2)])

    # Calculate pairs
    total_pairs = len(filled_sequence) // 2
    added_pairs = missing_length // 2

    return grouped_sequence, filled_sequence, total_pairs, added_pairs


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
expected_sequence_length = 4096

# Run the refined function
grouped_sequence, filled_sequence, total_pairs, added_pairs = identify_and_fill_gaps_properly(
    original_sequence, expected_sequence_length
)

# Display results
print("Original Sequence Length:", len(original_sequence))
print("Expected Sequence Length:", expected_sequence_length)
print("Filled Sequence Length:", len(filled_sequence))
print("Total Pairs in Filled Sequence:", total_pairs)
print("Added Pairs:", added_pairs)
print("\nGrouped Sequence (First 200 pairs):\n", grouped_sequence)  # Show first 200 pairs

# Optionally save the grouped sequence to a file
with open("grouped_sequence_refined.txt", "w") as file:
    file.write(grouped_sequence)
