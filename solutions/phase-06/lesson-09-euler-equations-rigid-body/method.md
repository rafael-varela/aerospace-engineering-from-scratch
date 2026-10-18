# P06.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. m = 10 + s/5; |r|² = 2.5 m²; ΔI = m(2.5E − rrᵀ). 2. Add to diag(260.4, 208.3, 135.4); only I₁₂ appears because r₃ = 0. 3. `principal_axes(I)`; the columns of R. 4. Angles arccos|R₁ₖ|; identify the middle moment.

## 2.
(a) ω · (2) ⇒ dT/dt = 0; (Iω) · (2) ⇒ d(H²)/dt = 0; H_i = CᵀIω and dH_i/dt = Cᵀ(d(Iω)/dt + ω × Iω) = 0. (b) H² − I₂(2T) term by term; set to zero. (c) ω_i · H_i = ω · Iω = 2T, with H_i fixed. (d) Insert and use (6): the first equation gives −A₁λ = (I₂ − I₃)A₂A₃/I₁, the second A₂λ = (I₃ − I₁)A₃A₁/I₂, the third −mA₃λ = (I₁ − I₂)A₁A₂/I₃; evaluate (4) at τ = 0 (sn = 0, cn = dn = 1) and at τ = K (sn = 1, cn = 0, dn = k′) to get A₁ and A₃, A₂; then λ and m. (e) K(0) = π/2; λ with a → 0; K(m) ≈ ln(4/k′) and λ → σ as m → 1.

## 3.
(a) `kinetic_energy`, `np.linalg.norm(I @ w)`, `polhode_type`. (b) `elliptic_parameters`: period, m. (c) `simulate` over one period; H_i at t = 0; for each sample the body axis k in inertial components is row k of `dcm_from_quat(q)`; arccos of its dot product with Ĥ. (d) H is fixed in space; the body axis cones around it.

## 4.
(a) spin = 2π(2 + s/10), ε = 10^−(1 + s mod 4); `elliptic_parameters(T_HANDLE, [ε spin, spin, 0])`; `intermediate_growth_rate`. (b) `simulate(…, 3 * period, n_samples=6001, rtol=1e-11, atol=1e-12)`; `flip_times(t, w[:, 1])`; `elliptic_solution`; `invariant_drift`. (c) Same at 10⁻⁶. (d) 60 / flip interval.

## 5.
(a) `intermediate_growth_rate(I, 0.5)`, `elliptic_parameters`; half the flip interval. (b) ln(0.5/10⁻⁶)/σ against (1/σ) ln(4/k′). (c) `simulate(I, w0, 200, n_samples=20001)`, `flip_times`. (d) Each factor of 10 adds ln 10/σ to the first reversal (half of 2 ln 10/σ).

## 6.
(a) The torque or the H transformation in the wrong frame; a quaternion that is not renormalised or a wrong product order in q̇. (b) Logarithmic sensitivity of the phase to the disturbance. (c) 301 between 300 and 500: intermediate, but with σ tiny; compare σ with the spin rate. (d) H fixed, ω moves. (e) The attitude needs another integral. (f) Rigidity (fuel slosh): look for slowly growing nutation or a changing tumble period.
