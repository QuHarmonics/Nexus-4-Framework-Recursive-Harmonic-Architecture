def recursive_dimension_solver(a, b, max_dim=5):
    print(f"{'Dim':<5} {'a^n':<10} {'b^n':<10} {'Sum':<10} {'C (root)':<10}")
    print("-" * 50)
    for n in range(1, max_dim + 1):
        a_pow = a ** n
        b_pow = b ** n
        total = a_pow + b_pow
        c = total ** (1/n)
        print(f"{n:<5} {a_pow:<10.4f} {b_pow:<10.4f} {total:<10.4f} {c:<10.4f}")
