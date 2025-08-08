    # … after spike-handling …

    # ONLY now test for solution
    c1 = truth[0] or (not truth[1]) or truth[2]
    c2 = (not truth[0]) or truth[1] or truth[2]
    c3 = truth[0] or truth[1] or (not truth[2])

    if spike.any() and c1 and c2 and c3:
        print(f"Solved at t={t}: x1={int(truth[0])}, "
              f"x2={int(truth[1])}, x3={int(truth[2])}")
        break
