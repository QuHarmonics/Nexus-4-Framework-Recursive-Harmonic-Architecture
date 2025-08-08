for _ in range(iterations):
    # Dynamic perturbation based on divergence
    divergence = np.linalg.norm(current_guess - target_lattice)
    perturbation = 0.1 * np.sin(divergence)
    
    # Refine the guess lattice using dynamic perturbation
    current_guess = np.abs(np.sin(current_guess + perturbation))
    guess_lattices.append(current_guess)

    # Update divergence measures
    divergence_measures.append(divergence)

    # Early stopping condition
    if divergence < 1e-6:
        print(f"Converged at iteration {_+1}")
        break
