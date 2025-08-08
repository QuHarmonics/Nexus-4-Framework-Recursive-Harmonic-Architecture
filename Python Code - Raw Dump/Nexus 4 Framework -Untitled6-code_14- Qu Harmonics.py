def navier_stokes_update_full(Psi, delta_t, viscosity=0.1, iterations=20):
    """
    Performs one full Navier-Stokes velocity field update step on 2D vector field Psi.
    Psi shape: (N, N, 2) representing velocity components (u,v).
    """
    N = Psi.shape[0]
    u = Psi[..., 0]
    v = Psi[..., 1]

    # Compute nonlinear advection terms using central difference
    def ddx(f):
        return (np.roll(f, -1, axis=1) - np.roll(f, 1, axis=1)) / 2
    def ddy(f):
        return (np.roll(f, -1, axis=0) - np.roll(f, 1, axis=0)) / 2

    u_x = ddx(u)
    u_y = ddy(u)
    v_x = ddx(v)
    v_y = ddy(v)

    adv_u = u * u_x + v * u_y
    adv_v = u * v_x + v * v_y

    # Diffusion (viscous term)
    laplace_u = np.roll(u, 1, axis=0) + np.roll(u, -1, axis=0) + np.roll(u, 1, axis=1) + np.roll(u, -1, axis=1) - 4 * u
    laplace_v = np.roll(v, 1, axis=0) + np.roll(v, -1, axis=0) + np.roll(v, 1, axis=1) + np.roll(v, -1, axis=1) - 4 * v

    u_new = u + delta_t * (viscosity * laplace_u - adv_u)
    v_new = v + delta_t * (viscosity * laplace_v - adv_v)

    # Pressure projection step to enforce incompressibility
    div = ddx(u_new) + ddy(v_new)
    p = np.zeros_like(u_new)

    # Solve Poisson equation ∇²p = div using Jacobi iteration
    for _ in range(iterations):
        p = (np.roll(p, 1, axis=0) + np.roll(p, -1, axis=0) + np.roll(p, 1, axis=1) + np.roll(p, -1, axis=1) - div) / 4

    # Subtract pressure gradient from velocity
    u_proj = u_new - delta_t * ddx(p)
    v_proj = v_new - delta_t * ddy(p)

    Psi_proj = np.stack([u_proj, v_proj], axis=-1)
    return Psi_proj
