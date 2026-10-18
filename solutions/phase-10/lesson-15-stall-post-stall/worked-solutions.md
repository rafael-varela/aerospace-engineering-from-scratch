# P10.15  Worked solutions

Each solution states: assumptions, units, the derivation or method, the result, at least one independent check, and the common mistakes. Numerical results come from `solution.py` and a probe script that calls it, given to the figures shown. The section curves are fits to the authentic Ladson data (NASA TM-4074); the wing results are computations.

## 1. The section model

(a) At 10° with a₀ = 5.741 /rad: f = 1 → **cl = 1.002**, f = 0.5 → **0.730**, f = 0 → **0.251**. A half-separated section keeps **72.9 %** of its attached-flow lift: ((1 + √0.5)/2)² = 0.729.

(b) At 45°: **cl = 1.000, cd = 1.000**; the resultant is √2 = **1.414 = 2 sin 45°**, perpendicular to the plate (the model has no chordwise force).

(c) Equation (1) is a small-angle result whose fully separated branch still grows linearly, giving 1.50 at 60° and 2.26 at 90°, while a plate's lift must fall to zero at 90°; the plate model's normal force 2 sin α is bounded and points the right way.

**Common mistakes.** Degrees in `kirchhoff_cl`.

## 2. Derivations

(a) f = 1: ((1 + 1)/2)² = 1, so cl = a₀α = 2πα, thin-airfoil theory. f = 0: ((1 + 0)/2)² = ¼.

(b) Minimise S(p) = (Xp − y)ᵀ(Xp − y): ∇S = 2Xᵀ(Xp − y) = 0 ⇒ **XᵀXp = Xᵀy**. With independent errors of variance σ² in y, p = (XᵀX)⁻¹Xᵀy has covariance (XᵀX)⁻¹Xᵀ(σ²I)X(XᵀX)⁻¹ = **σ²(XᵀX)⁻¹**, with σ² estimated by Σr²/(n − m). For a nonlinear model, linearise r(p) ≈ r(p*) + J(p − p*) near the solution: the same algebra holds with X → J. It is approximate (the curvature of the model is ignored) and it assumes independent random errors, which a model-form error such as a step-like stall violates.

(c) At convergence δ_new = δ, i.e. δ = α_eff − cl_sec(α_eff)/a₀ at every station, so the linear solver's local lift a₀(α_eff − δ) equals cl_sec(α_eff). The linear solver enforces 2bΣA_n sin nθ = ½c·a₀(α_eff − δ), which is therefore ½c·cl_sec(α_eff): equation (4). a₀ only changes the path, not the fixed point.

(d) With the exact solution A* = JA* + b, the error e = A − A* obeys e ← (1 − ω)e + ωJe = [(1 − ω)I + ωJ]e. Its eigenvalues are 1 + ω(λ − 1); the iteration converges for every start iff all have modulus below 1. For real λ < 0: 1 + ω(λ − 1) < 1 always (ω > 0), and > −1 iff **ω < 2/(1 − λ)**; the most negative λ sets the limit. For the trainer's rectangle (16 terms, a₀ = 5.74) the code gives 0.062, and the lesson's experiment confirms convergence at 0.05 and failure at 0.1.

**Common mistakes.** Quoting the fit's standard errors as the model's accuracy.

## 3. A tilt-wing in transition

q = 245 Pa, S = 12 m².

| α | cl | cd | lift (N) | drag (N) | normal force (N) |
|---|---|---|---|---|---|
| 35° | 0.940 | 0.658 | 2763 | 1934 | 3373 |
| 45° | 1.000 | 1.000 | 2940 | 2940 | 4158 |
| 60° | 0.866 | 1.500 | 2546 | 4410 | 5092 |

(b) The lift supports the aircraft (the rotors' thrust can be reduced by that much); the drag, as large as the lift or larger, decelerates it and must be overcome by the thrust. Blowing the wing with the propellers' slipstream raises the local dynamic pressure and lowers the local angle of attack, so the wing works in nearly attached flow during conversion instead of as a stalled plate.

(c) Steady flow: during conversion α changes in seconds, and dynamic-stall effects (`P10.19`) and the slipstream make the forces very different from a steady plate's.

**Common mistakes.** Taking lift along the wing's normal.

## 4. Washing out the tapered wing

(a)

| washout | α at first stall | 2|y|/b | C_L there | e at 5° |
|---|---|---|---|---|
| 0° | 17.10° | 0.556 | 1.384 | 0.983 |
| 2° | 18.07° | 0.471 | 1.392 | 0.966 |
| 4° | 18.79° | 0.290 | 1.381 | 0.847 |
| 6° | 19.35° | 0.195 | 1.358 | 0.569 |

(b) **4°** brings the first stall inboard of 30 % of the semi-span (0.290). The cost in C_L at stall onset is small (1.384 → 1.381) because the higher α makes up for the reduced tip lift; the stall angle rises by 1.7°.

(c) Twist adds a load distribution that does not scale with α (positive inboard, negative at the tips at the zero-lift angle), so at a low C_L like that at 5° the load is far from elliptic and e is low; at high C_L the twist share is smaller. e of a twisted wing depends on C_L; the right measure is C_Di against C_L over the operating range (in particular at cruise C_L).

**Common mistakes.** Applying washout as nose-up twist (sign).

## 5. A higher Reynolds number

(a) Re = 8.86 × 10⁶ fit: section stall angle **19.13°**, cl_max **1.700**.

(b)

| planform | first stall α | 2|y|/b | C_L | largest C_L in sweep (α) |
|---|---|---|---|---|
| rectangular | 21.15° | 0.000 | 1.585 | 1.591 (21.5°) |
| tapered 0.45 | 21.56° | 0.556 | 1.652 | 1.658 (22.0°) |
| 4° washout | 21.39° | 0.000 | 1.512 | 1.523 (23.5°) |

(One angle in each sweep did not converge, past stall.)

(c) The stations of first stall did not change (they are set by the planform's load distribution, not by the section); every C_L rose by 0.27–0.30 with the section's cl_max, and the gentler section stall makes the sweep's maximum slightly exceed the critical-section value. The trainer's wing reaches Re = 4.5 × 10⁶ at its 55 m/s cruise and 2.5 × 10⁶ at approach, never 8.86 × 10⁶: the Re = 2 × 10⁶ fit is the right choice for its stall, which happens at low speed.

**Common mistakes.** Keeping the old stall angle in `stall_onset_alpha`.

## 6. What can this model say about hysteresis?

(a) C_L on the way down after a sweep up to 22°:

| terms | 16° | 14° | 12° |
|---|---|---|---|
| 8 | 0.717 | 0.634 | 0.950 |
| 16 | 0.810 | 0.696 | 0.645 |
| 32 | 0.828 | 0.731 | 0.661 |

(On the way up the three resolutions agree to the truncation error: 0.950 at 12°.)

(b) Below stall the section curve has a positive slope and the lifting-line equation has one solution, which every resolution approximates better and better. Past stall a negative slope allows many solutions (stalled patches of any width); with more terms the series can represent narrower stalled regions, so the set of solutions, and which one the continuation lands on, changes with the resolution. The 8-term solution reattaches by 12°, the 16- and 32-term ones do not.

(c) "The model shows that a stalled branch can exist, so some hysteresis is plausible, but its size and even its existence at 12° depend on the number of Fourier terms, so it predicts nothing quantitative. Hysteresis of the real wing has to come from a wind-tunnel test or a model with spanwise and unsteady flow."

(d) (1) Spanwise flow and three-dimensional separation (a 3-D viscous model or at least a vortex-lattice method with a stall model and spanwise coupling); (2) time dependence: an unsteady wake and dynamic-stall section behaviour (lag of separation and reattachment, `P10.19`), with measured 2-D hysteresis data for the section.

**Common mistakes.** Reporting a post-stall C_L curve from one resolution as "the" result.

## Check-yourself answers

1. cl = a₀α((1 + √f)/2)²; f = 1: a₀α (2πα for a₀ = 2π); f = 0: a₀α/4.
2. How quickly the separation point moves forward with α (the width of the transition from attached to separated). The NACA 0012 at Re = 2 × 10⁶ stalls abruptly from the leading edge, so the best fit is nearly a step (S = 0.28°).
3. 45°: cl = 1, cd = 1. 90°: cl = 0, cd = 2.
4. It treats the induced angle explicitly; the high Fourier harmonics amplify errors (eigenvalues of the iteration matrix far below −1), so ω must be below about 0.06. The decambering iteration solves the linear lifting line, including the induced angle, implicitly in each sweep and only iterates on the section's nonlinearity.
5. The tip vortices induce the largest downwash near the tips of a rectangular wing, so its outer sections work at lower effective α and the root sections reach cl_max first.
6. It lowers the tips' angle of attack, moving the first stall inboard; it lowers C_L at stall onset at a given planform (1.31 → 1.22 for the rectangle) or leaves it nearly unchanged when it corrects a tapered wing's tip stall.
7. Past cl_max the section curve has a negative slope, so the lifting-line equation admits several distributions with different stalled patches; the solution found depends on the start and on the resolution, so the model cannot predict a hysteresis loop's size.
