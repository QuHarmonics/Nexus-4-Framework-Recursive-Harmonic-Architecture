from scipy.optimize import curve_fit

def exp_decay(t, A, B, tau):
    return A + B * np.exp(-t / tau)

popt, _ = curve_fit(exp_decay, list(range(len(M_series))), M_series)
A, B, tau = popt
print(f"Fitted τ (relaxation time): {tau:.2f}")