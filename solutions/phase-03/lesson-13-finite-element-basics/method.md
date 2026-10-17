# P03.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) h = L/2; k_e = (EA/h)[[1, −1], [−1, 1]] for both; K = (EA/h)[[1, −1, 0], [−1, 2, −1], [0, −1, 1]]; f = [0, 0, P]; fix u₀ = 0; solve the 2 × 2 system for u₁, u₂. (b) σ_e = E(u_{e+1} − u_e)/h; r₀ = (Ku − f)₀. (c) f = qh/2 × [1, 2, 1] (interior node gets both halves); solve; compare with qL²/(2EA); nodal exactness of the uniform bar.

## 2.
(a) ∫(−(EAu')' − b)v dx = 0; parts: ∫EAu'v' dx − [EAu'v]₀ᴸ; v(0) = 0 kills the lower limit, EAu'(L) = P gives Pv(L). Essential: u(0) = 0; natural: the tip force. (b) k_ij = ∫EAN_i'N_j' dx = EA(h/2)∫(∓1/h)(∓1/h) dξ; m_ij = ρA(h/2)∫N_iN_j dξ with N = (1 ∓ ξ)/2. (c) vᵀKv = ∫EA(Σv_jN_j')² dx ≥ 0; zero iff v_h' = 0 iff v constant; after fixing node 0 a constant admissible v is zero. (d) ∂Π/∂u = Ku − f; the Hessian is K, positive definite on the free dofs; substitute Ku = f into Π.

## 3.
(a) Replace σ by E(u' − αΔT) in ∫σAv' dx = ∫bv dx + Pv(L); move the thermal term to the right; evaluate ∫EAαΔT N_i' dx on each element. (b) `assemble`, then add EAαΔT[−1, 1] to f[e:e+2] for every element; `solve_bar` with {0: 0, 4: 0}; σ = E u' − EαΔT. (c) Same with {0: 0}. (d) Fixed–fixed with any ΔT(x): u ≡ 0 is the exact solution only if the thermal load is self-equilibrating... work it out: the exact σ is constant (equilibrium with b = 0), equal to −Eα × (mean ΔT); the finite-element solution with the consistent thermal load reproduces it because the exact u is piecewise linear (in fact linear) and lies in the trial space.

## 4.
(a) P = m_tipΩ²L; spar load ρAΩ²L²/2; (5) at x = 0 and x = L; reaction = −(spar load + P). (b) `solve_rotor_blade(case, 8)`; `exact_rotating_bar` at the nodes, at h/2 and at 0; exact: nodal u, reaction; approximate: element stress. (c) σ₀ + (σ₀ − σ₁)(x_c0 − 0)/(x_c1 − x_c0); compare. (d) `convergence_study`, `observed_order`.

## 5.
(a) `assemble(..., rho)`; `axial_frequencies` with M and with `np.diag(M.sum(axis=1))`; divide by `exact_axial_frequencies`. (b) Sign of ratio − 1 and its trend; p from log₂ of successive error ratios. (c) Mode 3 has 2.5 half-waves on L; find the n with |ratio − 1| < 0.01; elements per half-wave = n/2.5. (d) Compare 254.6 Hz with kΩ up to 66.7 Hz; bending frequencies are an order of magnitude lower and are the ones that cross the harmonics.

## 6.
(a) Loop n; `assemble` with A(x) = A₀(1 − cx/L); `element_stress`[-1] against P/A(L). (b) In the last element EA_mid u' = P exactly (the discrete equation of a bar with one load), so σ_e = P/A_mid = σ_exact(x_mid); the distance to the tip is h/2, so the error is ≈ σ'(L)h/2. (c) σ(L) = P/A(L) from the natural boundary condition; or linear extrapolation of the last two element stresses, exact for a linear σ (not quite here). (d) Rank: bending (pinned ends carry no moment but the strut deflects laterally), pin-hole stress concentration, axial-only section behaviour; consequences: lateral deflection under load, local strain at the lug, a measured axial strain that varies around the circumference.
