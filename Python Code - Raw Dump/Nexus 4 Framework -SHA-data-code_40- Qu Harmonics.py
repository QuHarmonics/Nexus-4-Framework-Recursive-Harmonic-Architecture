transitions = []

for x in range(255):  # skip 255 to avoid overflow
    x_next = x + 1
    r = hash_mod_residue(x, modulus)
    r_next = hash_mod_residue(x_next, modulus)
    delta_r = (r_next - r) % modulus

    b1 = format(x, '08b')
    b2 = format(x_next, '08b')
    bit_diff = sum(c1 != c2 for c1, c2 in zip(b1, b2))

    transitions.append({
        'x': x,
        'x+1': x_next,
        'residue': r,
        'residue_next': r_next,
        'delta_residue': delta_r,
        'bit_flip': bit_diff
    })

df_trans = pd.DataFrame(transitions)
display(df_trans.head(10))
