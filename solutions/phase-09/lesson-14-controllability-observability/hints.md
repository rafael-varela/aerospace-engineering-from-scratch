# P09.14  Hints

One hint per problem. Read only the one you need.

1. W_c = [B, AB] and W_o = [C; CA]. With A = [[0, 1], [0, 0]], AB shifts the second entry of B into the first. Only λ = 0 needs testing because it is the only eigenvalue.

2. (a) Every reachable state is a combination of e^{At}B, whose series can be rewritten with A⁰, …, A^{n−1} only. (b) Multiply ẋ = Ax + Bu on the left by w. (c) Differentiate F(t) = e^{At}BBᵀe^{Aᵀt} and integrate from 0 to T. (d) Write any other admissible input as u + v and show that the cross term ∫uᵀv dt vanishes. (e) Integrate e^{At}B = (t/J, 1/J)ᵀ times its transpose.

3. Build W for an arbitrary β by copying `wheel_axes` with the angle as a parameter. `np.linalg.eigvalsh(W @ W.T)` gives the torque authority per direction. cos²β = 1/3 is the angle at which the nominal pyramid is isotropic.

4. Use `ctrb`, `numerical_rank` and the PBH functions on the CW matrix. Physically, look at the CW equation for ÿ: radial thrust does not appear in it, so some combination of x and ẏ cannot be changed by it. Find that combination and say what it means for the orbit. For sensing, ask what happens to a chaser displaced purely along-track.

5. Bisect on ε with `is_controllable`. The Gramian eigenvalue scales as ε² in the weak direction. For (b), think of what "unobservable" means for an estimator's covariance and what a small Gramian eigenvalue means instead. For (c), think torque authority, momentum capacity and wheel friction.

6. Use `min_energy_input` for each configuration and evaluate u on a time grid. The body torque the slew needs points roughly along J·Δθ: take the dot product of that direction with each wheel axis. (u(t) = Bᵀ(...) lies in the row space of W, and the minimum-norm solution of Wu = τ is Wᵀ(WWᵀ)⁻¹τ.)
