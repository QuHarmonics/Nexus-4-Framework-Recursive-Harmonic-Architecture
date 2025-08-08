def adaptive_harmonic_fold(Psi, delta, sti, alpha_0=0.35, kappa=1.5):
    # Calculate adaptive gain per spatial location
    gain = alpha_0 * (1 + kappa * (1 - sti))
    gain = np.clip(gain, 0, 1)  # Prevent excessive gain
    correction = -gain[..., None] * delta
    return Psi + correction

for t in range(time_steps):
    Psi_new = navier_stokes_update(Psi, delta_t)
    delta = Psi_new - Psi
    sti = symbolic_trust_index(delta, prev_delta)
    
    unstable_mask = sti < stability_threshold
    delta_corrected = np.where(unstable_mask[..., None], delta, 0)
    
    Psi_new = adaptive_harmonic_fold(Psi_new, delta_corrected, sti)
    
    prev_delta = delta
    Psi = Psi_new
    
    print(f"Step {t}: Avg STI = {np.mean(sti):.4f}")
