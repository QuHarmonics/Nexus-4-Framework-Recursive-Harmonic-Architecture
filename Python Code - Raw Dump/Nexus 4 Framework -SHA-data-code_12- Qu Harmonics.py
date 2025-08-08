    # Δ–Σ spike detection
    spike = np.abs(err) > tau
    if spike.any():
        idx = np.where(spike)[0][0]
        var, bit = divmod(idx, 3)

        print(f"\nTick {t}: Spike at slot {idx} → var {var+1}, bit {bit}")
        print("  Before toggle truth:", truth.astype(int))

        if bit in (0, 2):
            truth[var] = not truth[var]
            print(f"  Toggled x{var+1} → now", int(truth[var]))
        else:
            print("  Spike on spacer—no toggle")

        phi[idx] = 0.0
        print("  φ at idx cleared; current φ:", np.round(phi, 2))

    # Now show clause satisfaction status
    c1 = truth[0] or (not truth[1]) or truth[2]
    c2 = (not truth[0]) or truth[1] or truth[2]
    c3 = truth[0] or truth[1] or (not truth[2])
    print(f"Tick {t}: Clause status – c1={c1}, c2={c2}, c3={c3}")
