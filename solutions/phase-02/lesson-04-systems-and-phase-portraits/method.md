# P02.04  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. Classify by hand

τ = a + d, Δ = ad − bc with c = −s/10. Δ > 0 rules out a saddle, τ < 0 gives stability, and the sign of τ² − 4Δ picks node (positive) or spiral (negative). Eigenvalues from λ = τ/2 ± √(τ²/4 − Δ); if complex, ω_n = √Δ, ζ = −τ/(2ω_n), ω_d = Im λ, period 2π/ω_d, t_half = ln 2/|Re λ|. Nullclines: ẋ₁ = 0 ⇒ x₂ = (−a/b)x₁; ẋ₂ = 0 ⇒ x₂ = (−c/d)x₁. Sense of rotation: on the positive x₁ axis ẋ₂ = c x₁, so c < 0 means the trajectory moves to negative x₂, i.e. clockwise. Check with `classify`, `eigen_structure`, `modal_parameters`, `nullcline_slopes`.

## 2. Derivations

(a) det(A − λI) = (a − λ)(d − λ) − bc = λ² − (a + d)λ + (ad − bc). Vieta gives the sum and product of roots.
(b) Let **z**(t) = **v**e^{λt} with **v** = **p** + i**q**, λ = σ + iω. Then **ż** = A**z**. Write **z** = **u** + i**w**; since A is real, **u̇** + i**ẇ** = A**u** + iA**w**, so **u̇** = A**u** and **ẇ** = A**w**. Expand e^{λt}(**p** + i**q**) with Euler's formula to read off **u** and **w**. In the basis (**p**, **q**) the coordinates are e^{σt}(cos ωt, −sin ωt): a circle in that basis, an ellipse in the standard one, radius scaled by e^{σt}.
(c) dQ/dt = 2c x₁ẋ₁ − 2a(ẋ₁x₂ + x₁ẋ₂) − 2b x₂ẋ₂. Substitute, collect x₁², x₁x₂, x₂² terms; each coefficient contains the factor (a + d) or vanishes identically. Definiteness: Q's matrix is [[c, −a], [−a, −b]] with determinant −bc − a² = Δ when d = −a; positive determinant means definite (ellipses), negative means indefinite (hyperbolas).
(d) Row 1 (speed): u̇ = −gθ − (2g/(L/D))(u/V₀). Row 2 (flight-path angle): V₀θ̇ = (lift excess)/m = g·2u/V₀. Divide row 1 by V₀ to use the state (u/V₀, θ). Then Δ = 2g²/V₀² and τ = −2g/(V₀ L/D); ω_n = √Δ, ζ = −τ/(2√Δ).

## 3. Two tanks

(a) τ = −(k₁ + k₂), Δ = k₁k₂ − k₂k₁ = 0: non-isolated. Every state with k₁m₁ = k₂m₂ is an equilibrium (no net transfer), so there is a whole line of them, one for each total mass.
(b) λ = 0 with eigenvector (k₂, k₁) (the line of equilibria); λ = −(k₁ + k₂) with eigenvector (1, −1) (mass moves from one tank to the other along lines of constant total).
(c) m(t) = c₁(k₂, k₁) + c₂(1, −1)e^{−(k₁+k₂)t}; fix c₁, c₂ from m(0) = (100, 0). Final masses are the λ = 0 component; the time constant is 1/(k₁ + k₂). Check with `solve_linear`.
(d) Total mass m₁ + m₂: (1, 1) is a left eigenvector of A with eigenvalue 0 (the column sums of A are zero), which is the same statement as Δ = 0 with a conserved linear functional.

## 4. Your pitch damping

Build the matrix with `short_period_matrix(m_q=...)`; `classify` and `modal_parameters` give (a). For (b), `solve_linear` on t = linspace(0, 8, 80001) with x₀ = (5°, 0) in radians; the peak pitch rate is the extreme of x[:, 1] (negative: the nose drops to shed the extra angle of attack); the settling time is the time just after the last sample with |α| ≥ 0.5°. (c) `rk4_trajectory(A, x0, 8, 800)` against `solve_linear` on the same grid, max absolute difference. (d) `phase_portrait` with starts on a circle; nullclines from `nullcline_slopes`; trajectories cross the α-nullcline vertically (α̇ = 0) and the q-nullcline horizontally (q̇ = 0).

## 5. Judgment

(a) Time-scale separation: the fast mode completes several cycles before the slow mode's state has changed appreciably, so holding the slow state constant is accurate; the approximations fail when the periods are within a factor of a few (sailplanes, human-powered aircraft, very light short periods). (b) Doubling in 1.6 s is slow compared with a pilot's 0.3–0.5 s reaction and very slow compared with a flight-control computer; unstable-but-flyable means the loop (pilot or computer) supplies the missing M_α. Certification for transports requires stability without augmentation only in specified failure cases, so airliners are designed stable; relaxed-stability designs rely on redundancy. (c) Candidates: M_α̇ (downwash lag) adds damping; fuselage flexibility can reduce the effective tail moment; compressibility changes M_q and M_α; the frozen-speed assumption itself. Say for each whether it raises or lowers ζ. (d) Any amplitude dependence of frequency or decay, a response asymmetric between nose-up and nose-down, or saturation of α.

## 6. Walking into the saddle

(a) τ = −2 is fixed. Spiral→node: τ² = 4Δ ⇒ Δ = 1 ⇒ M_α = 0.96 − 1. Node→saddle: Δ = 0 ⇒ M_α = 0.96. Sweep M_α with `classify` on a fine grid and find where the string changes.
(b) For each M_α: λ₊ = −1 + √(1 − Δ) with Δ = 0.96 − M_α; doubling time ln 2/λ₊; unstable eigenvector (b, λ₊ − a) = (1, λ₊ + 0.8), angle arctan(λ₊ + 0.8). As λ₊ grows the eigenvector tilts toward the q axis: a faster divergence needs more pitch rate per unit angle of attack.
(c) With m = M_q + k: τ = −0.8 + m, Δ = −0.8m − 2 (positive only for m < −2.5). ζ = 0.5 means τ² = Δ, a quadratic in m; compute its discriminant. Then minimise ζ² = τ²/(4Δ) over m < −2.5 by differentiation: the minimum is above 1, so every stabilised configuration is a node. The reason: adding pitch damping raises |τ| and Δ together, but Δ only through the small product −0.8m, so the parabola τ² = 4Δ is never reached from below.
(d) Set Δ = τ²/(4ζ²) with τ = −2 and ζ = 0.45, then M_α = 0.96 − Δ. Verify with `modal_parameters`. The conclusion concerns what the feedback must measure.
