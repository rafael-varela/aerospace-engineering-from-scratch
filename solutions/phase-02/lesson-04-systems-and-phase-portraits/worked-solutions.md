# P02.04  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem); the "independent check" line says how it was verified by a different route. Units: s⁻¹ for rates, s for times, degrees where stated; angles inside the code are radians.

## 1. Classify by hand

**Assumptions.** Constant-coefficient linear system; the state is small enough that no nonlinearity matters.

**Derivation.** A = [[−1, 2], [−1.5, −3]]. τ = −4, Δ = (−1)(−3) − (2)(−1.5) = 3 + 3 = 6. Δ > 0: not a saddle; τ < 0: stable; τ² − 4Δ = 16 − 24 = −8 < 0: spiral. **Stable spiral.** Eigenvalues λ = −2 ± √(4 − 6) = **−2 ± 1.414i s⁻¹**. ω_n = √6 = **2.449 rad/s**, ζ = 4/(2·2.449) = **0.816**, ω_d = √2 = 1.414 rad/s, period 2π/1.414 = **4.44 s**, t_half = ln 2/2 = **0.347 s**. Nullclines: ẋ₁ = 0 on x₂ = x₁/2 (slope **0.5**); ẋ₂ = 0 on x₂ = −x₁/2 (slope **−0.5**). On the positive x₁ axis ẋ₂ = −1.5x₁ < 0: the trajectory heads to negative x₂, so rotation is **clockwise**.

**Result.** (τ, Δ) = (−4, 6); stable spiral; λ = −2 ± 1.414i; ω_n 2.449, ζ 0.816, T 4.44 s, t_half 0.347 s; slopes 0.5 and −0.5; clockwise.

**Independent check.** `classify(A)` → 'stable spiral'; `eigen_structure` → −2 ± 1.41421i; `modal_parameters` → {omega_n 2.4495, zeta 0.8165, omega_d 1.4142, period 4.4429, t_half 0.3466}; `nullcline_slopes` → (0.5, −0.5); `solve_linear` from (1, 0) gives x₂ = −0.068 at t = 0.05 (clockwise).

**Common mistakes.** Using τ² − 4Δ before checking the sign of Δ (a saddle has τ² − 4Δ > 0 too); quoting the period as 2π/ω_n (4.44 s is 2π/ω_d; 2π/ω_n = 2.57 s); calling ζ = 0.816 "overdamped".

## 2. Derivations

**(a)** det(A − λI) = (a − λ)(d − λ) − bc = λ² − (a + d)λ + (ad − bc) = λ² − τλ + Δ. Roots λ₁, λ₂ satisfy (λ − λ₁)(λ − λ₂) = λ² − (λ₁ + λ₂)λ + λ₁λ₂, so λ₁ + λ₂ = τ and λ₁λ₂ = Δ. Consequences: Δ < 0 ⇔ real roots of opposite sign (saddle); Δ > 0 and τ < 0 ⇔ both stable.

**(b)** **z** = (**p** + i**q**)e^{σt}(cos ωt + i sin ωt) = e^{σt}[(**p** cos ωt − **q** sin ωt) + i(**p** sin ωt + **q** cos ωt)]. Since **ż** = A**z** with A real, taking real and imaginary parts gives **u̇** = A**u** and **ẇ** = A**w** with **u** = Re **z**, **w** = Im **z**. In the basis (**p**, **q**), **u** has coordinates e^{σt}(cos ωt, −sin ωt): a point moving on a circle of shrinking (σ < 0) or growing (σ > 0) radius; the linear map from that basis to the standard one turns the circle into an ellipse. **u** and **w** are independent because **p** and **q** are (if **q** = α**p** then **v** would be a real eigenvector of a complex eigenvalue, impossible for real A).

**(c)** With d = −a: Q̇ = 2c x₁ẋ₁ − 2a(ẋ₁x₂ + x₁ẋ₂) − 2b x₂ẋ₂. Substitute ẋ₁ = a x₁ + b x₂, ẋ₂ = c x₁ − a x₂:
x₁²: 2ca − 2ac = 0; x₁x₂: 2cb − 2a(b − a) − 2a(a) ... collect carefully: from 2c x₁ẋ₁: 2ac x₁² + 2bc x₁x₂; from −2a ẋ₁x₂: −2a² x₁x₂ − 2ab x₂²; from −2a x₁ẋ₂: −2ac x₁² + 2a² x₁x₂; from −2b x₂ẋ₂: −2bc x₁x₂ + 2ab x₂². Everything cancels: Q̇ = 0. (For general d the surviving term is −2b(a + d)x₂² plus −2a(a + d)x₁x₂ ... the general statement is Q̇ = −2(a + d)(a x₁x₂ + b x₂²), zero when τ = 0.) Q = **x**ᵀ[[c, −a], [−a, −b]]**x**; the matrix has determinant −bc − a² = Δ (since Δ = −a² − bc when d = −a). Positive determinant ⇒ definite ⇒ level curves are ellipses (center); negative ⇒ indefinite ⇒ hyperbolas (saddle).

**(d)** Flight-path angle: the excess lift mg·2u/V₀ supplies centripetal acceleration V₀θ̇, so θ̇ = (2g/V₀)(u/V₀). Speed: u̇ = −gθ − (2g/(L/D))(u/V₀) (gravity component along the path plus the excess drag 2D u/V₀ with D = mg/(L/D)). Dividing the first... writing the state as (u/V₀, θ): d(u/V₀)/dt = −(2g/(V₀ L/D))(u/V₀) − (g/V₀)θ; dθ/dt = (2g/V₀)(u/V₀). So A = [[−2g/(V₀ L/D), −g/V₀], [2g/V₀, 0]], Δ = 2g²/V₀², τ = −2g/(V₀ L/D), ω_n = √2 g/V₀, ζ = −τ/(2ω_n) = (2g/(V₀ L/D))/(2√2 g/V₀) = 1/(√2 L/D). For V₀ = 230, L/D = 17: ω_n = 0.0603 rad/s, period 104.3 s, ζ = 0.0416, as in `main()`.

**Common mistakes.** Dropping the factor 2 in the lift excess (lift ∝ V²); treating ζ as depending on speed (it does not in this model); claiming **u** and **w** might be dependent.

## 3. Two tanks

**Assumptions.** Linear transfer rates; no losses; k₂ = 0.005 + 15/2000 = **0.0125 s⁻¹**, k₁ = 0.010 s⁻¹.

**Derivation.** A = [[−0.010, 0.0125], [0.010, −0.0125]]. τ = −0.0225, Δ = 0.010·0.0125 − 0.0125·0.010 = **0**: **non-isolated** fixed point. Every state with 0.010 m₁ = 0.0125 m₂ (m₁/m₂ = 1.25) is an equilibrium: the flows balance, and there is one such state for each total mass, forming a line through the origin. Eigenvalues 0 and τ = **−0.0225 s⁻¹**. For λ = 0: (A)**v** = 0 ⇒ **v** ∝ (k₂, k₁) = (0.0125, 0.010) ∝ (0.781, 0.625): the line of equilibria. For λ = −0.0225: (A + 0.0225I)**v** = 0 ⇒ (0.0125, 0.0125; 0.010, 0.010)**v** = 0 ⇒ **v** ∝ (1, −1): the direction of approach (mass leaves one tank and enters the other, total fixed). Solution: m = c₁(k₂, k₁) + c₂(1, −1)e^{−0.0225t}; at t = 0: c₁k₂ + c₂ = 100, c₁k₁ − c₂ = 0 ⇒ c₁ = 100/(k₁ + k₂) = 4444.4, c₂ = 44.44. So **m₁(t) = 55.56 + 44.44 e^{−0.0225t}**, **m₂(t) = 44.44 − 44.44 e^{−0.0225t}** kg; final masses **55.56 and 44.44 kg**; time constant 1/(k₁ + k₂) = **44.4 s**. Conserved: m₁ + m₂ = 100 kg; (1, 1)A = 0, so (1, 1) is a left eigenvector with eigenvalue 0, the dual statement of Δ = 0.

**Result.** `solve_linear` at t = 0, 100, …, 600 s: m₁ = 100, 60.24, 56.05, 55.61, 55.56, 55.56, 55.56 kg; m₂ = 0, 39.76, 43.95, 44.39, 44.44, 44.44, 44.44 kg.

**Independent check.** 100e^{−2.25}·0.4444 + 55.56 = 60.24 at 100 s by hand; m₁ + m₂ = 100.00 at every tabulated time.

**Common mistakes.** Calling Δ = 0 "a center" or "unstable" (it is a stable line: trajectories reach it and stop); normalising the λ = 0 eigenvector as (k₁, k₂) instead of (k₂, k₁); quoting the time constant as 1/k₁ or 1/k₂.

## 4. Your pitch damping

**Assumptions.** Short-period approximation (speed frozen), small angles, M_q = −(0.6 + 15/25) = **−1.2 s⁻¹**, which for s = 15 is the lesson's own aircraft.

**Result.** `classify` → stable spiral; `modal_parameters` → ω_n **2.227 rad/s**, ζ **0.449**, period **3.157 s**, t_half **0.693 s**. Propagating (5°, 0) with `solve_linear` on 80 001 points over 8 s: peak pitch rate **−5.15°/s at t = 0.555 s** (nose-down, shedding the excess angle of attack); |α| < 0.5° for all t > **2.00 s**. `rk4_trajectory` with 800 steps differs from the closed form by at most **2.5 × 10⁻¹⁰ rad** (order-4 error at h = 0.01 s). Phase portrait: eight starts on a circle of radius 0.3 spiral clockwise into the origin in about two turns; the α-nullcline q = 0.8α is crossed vertically (α̇ = 0) and the q-nullcline q = −3.33α horizontally (q̇ = 0).

**Independent check.** The envelope e^{−t} drops 5° to 0.5° at t = ln 10 = 2.30 s; the actual crossing at 2.00 s is earlier because the oscillation is near a zero of cos(ω_d t − φ) there. Peak q from the closed form q(t) = −(Δ/ω_d) α₀ e^{−t} sin ω_d t with Δ/ω_d·α₀ = 4.96/1.99·5° = 12.46°/s has its extreme where tan ω_d t = ω_d ⇒ t = 0.555 s, q = −12.46·e^{−0.555}·sin(1.105) = −5.15°/s.

**Common mistakes.** Reporting the first time |α| drops below 0.5° (1.5 s) rather than the last crossing; peak q in rad/s; forgetting that q starts at zero so its first extreme is negative.

## 5. Judgment

**(a)** Periods 3.2 s and 104 s differ by a factor 33. During one short-period cycle the phugoid state (speed, flight-path angle) changes by 1/33 of its own cycle, i.e. a few percent, so freezing speed costs a few percent in the short-period eigenvalues; conversely the short period has decayed completely (t_half 0.7 s) long before the phugoid has moved, so the angle of attack is at its quasi-static value during the phugoid. The separation fails for very light, slow aircraft (sailplanes, human-powered aircraft, small UAVs) whose short period can be 10 s against a 30 s phugoid: then the 4×4 system is needed and the modes interact.

**(b)** The pilot feels the nose start to wander and then accelerate away in pitch, doubling every 1.6 s; he can counter it (reaction 0.3–0.5 s is well inside the doubling time) but must do so continuously, which is fatiguing and dangerous in turbulence or when distracted. "Unflyable" means a divergence faster than the loop closing it; a flight-control computer closes the loop at 50–100 Hz with effectively no delay and makes a 1.6 s doubling time a non-event, which is how relaxed-stability fighters fly. Airliners are nonetheless designed stable because certification requires safe flight after augmentation failures, and because the margin against the saddle must survive 20 % uncertainty in M_α.

**(c)** M_α̇ (lag of the downwash at the tail): a positive-rate α̇ produces a nose-down moment, adding damping; omitting it makes the model *under*-predict ζ, so it cannot explain a measured 0.35 against a predicted 0.45. Fuselage flexibility reduces the tail's effective incidence change and so the effective |M_q| and |M_α|, lowering ζ (the right direction). Compressibility at cruise Mach can reduce M_q relative to the low-speed estimate, lowering ζ. Unsteady aerodynamics (lift lag) also reduce effective damping. The frozen-speed assumption itself contributes at the few-percent level.

**(d)** The decay rate or frequency depends on the pulse size; the response to a nose-up pulse differs from a nose-down one; α saturates at a value (stall) regardless of pulse amplitude; a second frequency appears.

## 6. Walking into the saddle

**(a)** τ = −0.8 − 1.2 = −2 for all M_α; Δ = (−0.8)(−1.2) − M_α = 0.96 − M_α. Spiral → node at τ² = 4Δ: 4 = 3.84 − 4M_α ⇒ **M_α = −0.04 s⁻²**. Node → saddle at Δ = 0: **M_α = 0.96 s⁻²**. `classify` on steps of 0.25: stable spiral for M_α ≤ −0.25, stable node for 0 ≤ M_α ≤ 0.75, saddle for M_α ≥ 1.0, and on a 0.001 grid the strings change at −0.040 and 0.960.

**(b)** λ₊ = −1 + √(1 − Δ) = −1 + √(0.04 + M_α); unstable eigenvector (1, λ₊ + 0.8) → angle arctan(λ₊ + 0.8).

| M_α (s⁻²) | λ₊ (s⁻¹) | doubling time (s) | eigenvector angle (deg from α axis) |
|---|---|---|---|
| 1.0 | 0.0198 | 35.0 | 39.3 |
| 1.5 | 0.2410 | 2.88 | 46.1 |
| 2.0 | 0.4283 | 1.62 | 50.8 |
| 3.0 | 0.7436 | 0.93 | 57.1 |

The angle increases with λ₊ because along the unstable eigenvector q = (λ₊ + 0.8)α: a faster divergence needs more pitch rate per unit angle of attack (α̇ = λ₊α = −0.8α + q).

**(c)** Let m = M_q + k. τ = −0.8 + m, Δ = (−0.8)m − 2. A spiral needs Δ > 0, i.e. m < −2.5. ζ = 0.5 requires τ² = 4ζ²Δ = Δ: (m − 0.8)² = −0.8m − 2 ⇒ m² − 0.8m + 2.64 = 0, discriminant 0.64 − 10.56 = **−9.92 < 0**: no real m. More: ζ² = (m − 0.8)²/(4(−0.8m − 2)); with u = −m > 2.5, d/du of (u + 0.8)²/(0.8u − 2) vanishes at 0.8u = 4.64, u = 5.8, where ζ² = 6.6²/(4·2.64) = 4.125, **ζ_min = 2.03** at M_q + k = −5.8. Every stabilised configuration is an overdamped node: pitch damping alone cannot make a statically unstable aircraft oscillate nicely, because feedback on q raises |τ| linearly while raising Δ only through the weak product (Z_α/V₀)M_q; the point (τ, Δ) moves along a line that never crosses into the spiral region.

**(d)** With τ = −2 and ζ = 0.45: Δ = τ²/(4ζ²) = 4/0.81 = 4.938, so M_α,eff = 0.96 − 4.938 = **−3.98 s⁻²**. `modal_parameters(short_period_matrix(m_alpha=-3.978))` → ζ = 0.450, ω_n = 2.222 rad/s. The augmentation system must supply the whole static stability itself, i.e. it must *sense angle of attack* (or an equivalent such as normal acceleration) and feed it to the elevator; a pitch-rate gyro alone is not enough. The lesson's "acceptable" short period at M_α = −4 is therefore also the design target of an α-feedback system for a relaxed-stability airframe.

**Independent check.** Each (τ, Δ) pair above substituted into λ² − τλ + Δ = 0 reproduces the quoted eigenvalues; `classify(short_period_matrix(m_alpha=2.0, m_q=-5.8))` → 'stable node' with eigenvalues −0.5 and −6.1, ratio consistent with ζ = 2.03 via (λ₁ + λ₂)²/(4λ₁λ₂) = 4.125.

**Common mistakes.** Setting Δ = 1 for the node boundary but forgetting τ² = 4Δ needs τ at the same M_q; concluding that "k < 0 enough" gives any ζ; in (d) forgetting that the required M_α is the *total* effective derivative, not the increment.

## Check-yourself answers (lesson.md)

1. λ² − τλ + Δ; Δ < 0 means a saddle (real eigenvalues of opposite sign), Δ > 0 means both eigenvalues have the same sign of real part.
2. τ = −3, Δ = −4 + 6 = 2 > 0, τ² − 4Δ = 9 − 8 = 1 > 0: stable node (eigenvalues −1 and −2).
3. ω_n = √(0.25 + 9) = 3.041 rad/s; ζ = 0.5/3.041 = 0.164; period 2π/3 = 2.09 s; t_half = ln 2/0.5 = 1.39 s.
4. Cayley–Hamilton gives (A − λI)² = 0, so the series for e^{At} truncates: **x**(t) = e^{λt}[I + (A − λI)t]**x**₀, and the t survives because there is only one eigenvector to carry a pure exponential.
5. The slow eigenvector (smaller |λ|): the fast component dies first, leaving the slow one.
6. Q = c x₁² − 2a x₁x₂ − b x₂²; for the oscillator it is −(k/m)x² − v² = −(2/m)·(mechanical energy).
7. M_α = +0.96 s⁻² (Δ = 0); beyond it a nose-up disturbance grows exponentially along the unstable eigenvector: static instability, the centre of gravity behind the neutral point.
8. ω_n = √2·9.81/100 = 0.139 rad/s, period 45.3 s; ζ = 1/(√2 L/D) depends only on the lift-to-drag ratio.
