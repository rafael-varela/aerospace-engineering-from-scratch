# P02.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. ω_n = √(k/m); f = ω_n/2π. 2. ζ = c/(2√(km)); ω_d = ω_n√(1 − ζ²); roots −ζω_n ± iω_d; regime from ζ. 3. x(t) = e^{−ζω_n t}[x₀ cos ω_d t + (ζω_n x₀/ω_d) sin ω_d t] with x₀ = 0.005 m; evaluate. 4. τ = 1/(ζω_n); δ = 2πζ/√(1 − ζ²). 5. c_crit = 2√(km). 6. `free_response(np.array([0.1, 0.5]), m, c, k, 0.005, 0.0)`.

## 2.
(a) Substitute, factor out e^{rt}, divide by m, complete the square or use the quadratic formula with c = 2ζω_n m, k = ω_n² m. (b) e^{(−σ ± iω_d)t} = e^{−σt}(cos ω_d t ± i sin ω_d t); real and imaginary parts are solutions; x(0) gives A = x₀; ẋ(0) = −σA + ω_d B gives B. (c) ẋ = (u̇ − ω_n u)e, ẍ = (ü − 2ω_n u̇ + ω_n²u)e; sum with coefficients 1, 2ω_n, ω_n²; all u and u̇ terms cancel. (d) Ẇ = x₁ẍ₂ − ẍ₁x₂; substitute ẍ_i = −(c/m)ẋ_i − (k/m)x_i. (e) Write x₀ cos + b sin as R cos(ω_d t − φ); peaks where cos = 1, spaced 2π/ω_d; ratio e^{σ·2π/ω_d}; take the log; solve the resulting equation for ζ.

## 3.
(a) m → L, c → R, k → 1/C: ω_n = 1/√(LC), ζ = (R/2)√(C/L). (b) ζ against 1; R_crit = 2√(L/C). (c) Underdamped: ringing at ω_d inside an e^{−Rt/2L} envelope; critical: single smooth decay. Filter designers usually want no ringing but a fast settle: ζ near 0.7 (Butterworth, met again in P02.08). (d) `natural_frequency(L, 1/C)`, `damping_ratio(L, R, 1/C)`, `classify_damping(L, R, 1/C)`, `damper_for_ratio(L, 1/C, 1.0)`.

## 4.
(a) Loop over the six ζ; `t = np.linspace(0, 3, 3001)`; `x, v = touchdown_response(t, m, k, z, v0, g)`; `max_stroke(m, k, z, v0, g)`; `strut_force(...).max()`; `loses_contact(x)`; `settling_time(t, x - xs, 0.02*xs)`. (b) `minimize_scalar(lambda z: strut_force(t, m, k, z, v0, g).max(), bounds=(0.05, 1), method="bounded")`; sweep `np.arange(0.1, 0.6, 0.001)` and take the first ζ with no bounce. (c) Pick a ζ just above the bounce threshold and near the force minimum; quote stroke, peak force and settling time. (d) `solve_ivp(lambda t, y: to_first_order(m, c, k) @ y, (0, 3), [-xs, v0], rtol=1e-10, atol=1e-12, t_eval=t)` against `free_response`.

## 5.
(a) Orifice damping: F_d = c₂ ẋ|ẋ| → 0 as ẋ → 0 at both the touchdown instant and maximum stroke; the force peak moves to mid-stroke; the "optimum ζ" of the linear model is an artefact of the spike. (b) The strut extends to its stop; the tyre, which can only push, unloads; the aircraft is airborne again with whatever vertical velocity it had. The tyre spring (stiffer, lightly damped) dominates the rebound. (c) Stroke: the gas-spring stiffening (it bounds the stroke). Load: the orifice damper (it sets the force history). (d) A force record with F(0) ≈ 0 and a peak near maximum stroke, or a force-velocity plot that is a parabola rather than a line.

## 6.
(a) ẋ = e^{−σt}[(−σx₀ + bω_d) cos ω_d t + (−σb − x₀ω_d) sin ω_d t]; set to zero; divide by cos; if the arctan comes out negative add π/ω_d. Evaluate with ω_n = 11.18, x₀ = −0.0785, v₀ = 3. (b) m Ẍ = mg − kX − cẊ = mg − F ⟹ F = m(g − Ẍ). Maximum F ⟺ minimum Ẍ (largest deceleration). For ζ = 0.3 the deceleration peaks shortly before the stroke peaks (sin and cos components); for ζ = 1.0 at t = 0 (the damper spike). Use `np.gradient` on the velocity to locate it. (c) F(0) = cv₀ = 2ζ√(km) v₀ grows with ζ and v₀; the force near maximum stroke ≈ k x_max falls with ζ. Halving v₀ halves the spike but reduces x_max less than proportionally, so the crossing moves to higher ζ: predict "up", then rerun `minimize_scalar` with v₀ = 1.5.
