# P02.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. Your own chute ($s = 15$: x_f = 4.5 m, y_f = 1.5 m)

**Derivation.** (a) r = 3. f(π) = π/2 − 3 < 0, so θ_f > π. Bisection on (π, 2π): f(4) = (4 − (−0.757))/(1 − (−0.654)) − 3 = 4.757/1.654 − 3 = −0.124; f(4.1) = (4.1 + 0.818)/(1 + 0.575) − 3 = +0.122; f(4.05) = (4.05 + 0.789)/(1 + 0.615) − 3 = −0.004: **θ_f ≈ 4.05 rad**. `cycloid_parameter(4.5, 1.5)` gives θ_f = 4.0516 rad (232.1°) and **k = 0.9295 m**.

(b) T = 4.0516 × √(0.9295/9.80665) = **1.2474 s**. Ramp: L = √(4.5² + 1.5²) = 4.743 m, T = 4.743 × √(2/(9.80665 × 1.5)) = **1.7490 s**. Saving: 1 − 1.2474/1.7490 = **28.7 %** (the ramp is 40.2 % slower).

(c) θ_f > π: yes. Lowest point y = 2k = **1.859 m** (0.36 m below the destination) at x = πk = **2.920 m**.

(d) `travel_time` on 401 points: 1.247388 s against 1.247387 s from (7): discretisation error **1.3 × 10⁻⁶ s** (chords are slightly slower than arcs; the error falls four-fold per doubling of points).

**Check.** The lesson's chute (r = 2) had θ_f = 3.508; a flatter geometry (r = 3) must have a larger θ_f ✓, and a longer flat run should widen the gap to the ramp: 40 % against 25 % ✓.

**Common mistakes.** Bracketing in (0, π) when f(π) < 0; using y_f = 1.5 in the ramp formula's length but forgetting it in the acceleration.

## 2. From the functional to the cycloid, and Snell

(a) J[y + εη] = ∫ F(x, y + εη, y′ + εη′) dx. d/dε at ε = 0: δJ = ∫ (F_y η + F_{y′} η′) dx. Integrate the second term by parts: ∫ F_{y′} η′ dx = [F_{y′} η]ₐᵇ − ∫ (d/dx F_{y′}) η dx. **The end conditions η(a) = η(b) = 0 kill the bracket.** So δJ = ∫ (F_y − d/dx F_{y′}) η dx, and since this vanishes for every such η, the fundamental lemma gives F_y − d/dx F_{y′} = 0.

(b) d/dx (F − y′F_{y′}) = (F_x + F_y y′ + F_{y′} y″) − y″F_{y′} − y′ d/dx(F_{y′}) = F_x + y′[F_y − d/dx F_{y′}]. On an extremal the bracket vanishes, leaving F_x. If F has no explicit x, F_x = 0 and F − y′F_{y′} is constant. If F_x ≠ 0 the quantity drifts at the rate F_x and is not conserved (as energy is not conserved for a time-dependent Lagrangian).

(c) F = (1 + y′²)^{1/2}(2gy)^{−1/2}, F_{y′} = y′(1 + y′²)^{−1/2}(2gy)^{−1/2}. F − y′F_{y′} = (2gy)^{−1/2}(1 + y′²)^{−1/2}[(1 + y′²) − y′²] = [2gy(1 + y′²)]^{−1/2} = C. Squaring and inverting: y(1 + y′²) = 1/(2gC²) ≡ 2k. Then y′² = (2k − y)/y. With y = k(1 − cos θ) = 2k sin²(θ/2): 2k − y = 2k cos²(θ/2), so y′ = cot(θ/2). dy = k sin θ dθ, and dx = dy/y′ = k sin θ tan(θ/2) dθ = 2k sin(θ/2)cos(θ/2) tan(θ/2) dθ = 2k sin²(θ/2) dθ = k(1 − cos θ) dθ. Integrating from θ = 0 with x = 0: x = k(θ − sin θ).

(d) ds² = dx² + dy² = k²[(1 − cos θ)² + sin²θ] dθ² = 2k²(1 − cos θ) dθ² = 4k² sin²(θ/2) dθ², so ds = 2k sin(θ/2) dθ and sin θ_path = dx/ds = k(1 − cos θ)/(2k sin(θ/2)) = sin(θ/2). The speed is v = √(2gy) = √(4gk sin²(θ/2)) = 2√(gk) sin(θ/2). Hence sin θ_path/v = 1/(2√(gk)) = 1/√(4gk), constant: Snell's law sin θ/c = const for a medium whose "speed of light" is √(2gy).

**Common mistakes.** Dropping the boundary term without saying why it vanishes; in (c), losing the factor 2 between k and the constant (any constant works, but then the cycloid has radius C/2).

## 3. A sound ray in a temperature inversion ($s = 15$: a = 0.015 s⁻¹, c₀ = 340 m/s, D = 2000 m)

**Assumptions.** Still air, sound speed a function of height only, geometric acoustics (wavelength ≪ 22 m, the arc height).

**Derivation.** (a) Snell: sin θ/c is constant. At the apex the ray is horizontal, θ = 90°, c = c₀ + ah, so the constant is 1/(c₀ + ah); at the ground c = c₀, giving sin θ₀ = c₀/(c₀ + ah). For a circle of radius R centred at height −c₀/a, the tangent at height y makes angle θ with the vertical where sin θ = (y + c₀/a)/R = (c₀ + ay)/(aR) = c(y)/(aR): sin θ/c = 1/(aR) is constant along the circle, so circles centred on y = −c₀/a are the extremals; the apex condition fixes R = (c₀ + ah)/a.

(b) Source and receiver on the ground at x = 0 and x = D: by symmetry the centre is at (D/2, −c₀/a) with c₀/a = 22 667 m, so R = √(1000² + 22 667²) = **22 689 m**; h = R − c₀/a = **22.05 m**; sin θ₀ = c₀/(c₀ + ah) = 340/340.331, θ₀ = 87.47°, launch elevation **2.53°** above the horizontal.

(c) `path_time` on the 2001-point arc: **5.88045 s**; on the straight ground path: **5.88235 s** (= 2000/340 exactly). The arc is faster by **1907 μs** although it is longer (by 0.65 m), because it spends its length where c is up to 0.33 m/s larger. The difference is small because the speed contrast is small (0.1 %), and first-order in the speed contrast: the arc saves time ∝ (ah/c₀) × D/c₀.

(d) `snell_invariant` on the arc has a spread of 10⁻¹⁵ around 1/(aR) = 2.9383 × 10⁻³ s/m. The ground path has sin θ = 1 and c = c₀ everywhere, so sin θ/c is also constant along it, but it is not an extremal: the Euler–Lagrange equation for F = √(1 + y′²)/c(y) contains F_y = −c′(y)√(1 + y′²)/c² = −a/c² ≠ 0 at y = 0, while d/dx F_{y′} = 0 for a straight path, so (2) is violated. Physically, a ray starting horizontally in a medium whose speed increases upward curves downward immediately and would go underground: the ground path is a boundary of the admissible region, and at a boundary the first variation need not vanish.

**Check.** Dimensions: c₀/a = m/s ÷ 1/s = m ✓. a → 0: R → ∞, h → 0, the arc becomes the straight line and the time difference → 0 ✓. (This is why aircraft noise carries further on a cold still night: rays that would leave upward are bent back to the ground.)

**Common mistakes.** Measuring θ from the horizontal; centring the circle at the ground.

## 4. Ritz on your chute ($s = 15$: x_f = 4.5 m, y_f = 1.5 m; cycloid 1.24739 s)

| modes | T [s] | excess |
|---|---|---|
| 1 | 1.29952 | 4.18 % |
| 2 | 1.28356 | 2.90 % |
| 4 | 1.26572 | 1.47 % |
| 8 | 1.25661 | 0.74 % |

(b) Monotone, and it must be: the family with n + 1 modes contains the n-mode family (set a_{n+1} = 0), so its minimum cannot be higher. Every entry is above 1.24739 s because the minimum over a subset of admissible paths is at least the minimum over all of them.

(c) The 4-mode Ritz path has slope 1.89 at x = 0; the cycloid's slope there is infinite (y ∝ x^{2/3}). A finite sine sum has a finite slope at the origin, so it can approach the vertical start only by piling up high modes, and each mode buys little: the excess roughly halves per *doubling* of modes.

(d) Family y = y_f x/x_f + a₁ x^{2/3}(1 − x/x_f) + a₂ sin(πx/x_f): Nelder–Mead gives a₁ = 1.409, a₂ = −0.157 and T = 1.24807 s, excess **0.055 %**: two parameters with the right singularity beat eight without it by a factor of 13.

**Check.** The 1-mode result is close to the lesson's 3.9 % for a similar geometry ✓.

**Common mistakes.** Starting Nelder–Mead far from the ramp (it can wander into paths that rise above the start, which the penalty rejects); comparing excesses to the ramp time instead of the cycloid time.

## 5. Stationary or not? ($s = 15$)

(a) Cycloid (x_f = 4.5, y_f = 1.5): ΔT = 9.66 × 10⁻⁴, 2.46 × 10⁻⁴, 6.22 × 10⁻⁵ s for ε = 0.1, 0.05, 0.025 m: ratios 3.92, 3.96, **exponent 1.97**. Circular arc of radius (4.5² + 1.5²)/9 = 2.5 m (time 1.2954 s): ΔT = 9.11 × 10⁻³, 4.45 × 10⁻³, 2.20 × 10⁻³ s: ratios 2.05, 2.03, **exponent 1.03**.

(b) J[y + εη] − J[y] = ε δJ + ½ε² δ²J + O(ε³). On the cycloid δJ = 0 (it is an extremal), so the change is quadratic: exponent 2, and the positive sign says δ²J > 0 in this direction, a minimum. On the arc δJ ≠ 0 and the linear term dominates: exponent 1, and a perturbation of the right sign *lowers* the time (ΔT > 0 here because η = sin(πx/x_f) lifts the arc, which is already too deep; −ε would lower T).

(c) x_f = 3.5 m, y_f = 2 m: x_f/y_f = 1.75 > π/2 = 1.571, so θ_f > π and the track dips. `cycloid_parameter(3.5, 2.0)`: θ_f = 3.309 rad, k = 1.0070 m; lowest point **2.014 m** deep at x = 3.164 m, so the part climbs 1.4 cm over the last 0.34 m.

(d) Time is length divided by speed, and speed is bought with depth: a track that dives deeper in the middle runs the long flat part of the journey faster, and the small climb at the end costs less time than the extra speed saved along the way. Nothing forbids an extremal from converting kinetic energy back to height; the functional only counts time.

**Common mistakes.** Taking |ΔT| without noticing the sign; concluding that exponent 2 alone proves a minimum (it shows stationarity plus one positive second variation).

## 6. Building the chute

(a) **Friction** (largest): work ≈ μmg × horizontal run = 0.2 × 4 m = 0.8 m of head against y_f = 2 m, 40 % of the energy, cutting final speed by 23 % and lengthening both times, but the cycloid more than the ramp (its steep start carries a large normal force in the curved part and loses its vertical-start advantage): the saving falls from 25 % toward 10–15 %. **Point mass vs rolling ball**: rescales g by 5/7 (solid sphere), lengthening every time by √(7/5) = 1.18 but leaving every *ratio* unchanged: no effect on the saving. **Released from rest**: an initial speed v₀ = 1 m/s moves the virtual cusp 5 cm above the start and slightly flattens the optimum: a per cent. **Planar track**: none for a guided part.

(b) Start with a straight drop at the friction angle (steeper than arctan μ = 11.3° so the part moves at all; in practice 45–60°) for the first 0.2–0.3 m, then join the cycloid tangentially. The cost is a few per cent of time relative to the ideal, far less than the 25 % the straight ramp loses; the gain is that the part starts at all.

(c) Two tracks between the same end points (4 m along, 2 m down), one cycloidal, one straight, a steel ball and two light gates. Predicted time ratio 1.253 for both a sliding mass and a rolling ball (g rescales equally). Measure five runs each; a ratio below about 1.15 means friction or rolling resistance is different between the tracks (rougher bends, the ball slipping on the steep start) and the frictionless model is not describing the chute.

(d) "The Euler–Lagrange equation proves the cycloid is *stationary*: no small change of shape alters the time to first order. It does not say minimum: a saddle would satisfy it too. The lesson needed the positive second variation (ΔT ∝ +ε²), the competing paths (ramp 25 % slower, arc 3 % slower, twenty random paths all slower) and the Ritz upper bounds converging from above to say *fastest*."

## Check-yourself answers

1. F_y − d/dx(F_{y′}) = 0. Integration by parts moved the derivative off the variation η (so that η appears as a plain factor and the fundamental lemma applies) and produced a boundary term that the end conditions kill.
2. If ∫ G η dx = 0 for every smooth η vanishing at the ends, then G ≡ 0. Without it one could only say the *integral* of G η vanishes, not that the integrand does.
3. F − y′F_{y′} = constant on extremals when F has no explicit x dependence. With x as time and F as the Lagrangian it is (minus) the energy, the Hamiltonian.
4. y′ = (dy/dθ)/(dx/dθ) = sin θ/(1 − cos θ) = cot(θ/2); 1 + y′² = 1/sin²(θ/2); y(1 + y′²) = k(1 − cos θ)/sin²(θ/2) = 2k sin²(θ/2)/sin²(θ/2) = 2k ✓.
5. Speed is √(2gy): the only way to be fast early is to be deep early, and a vertical drop gains depth at the greatest rate per unit path length.
6. x_f/y_f = 3 > π/2, so θ_f > π: yes, it dips.
7. Quadratic: a quarter, 2.5 × 10⁻⁴ s. On the ramp the change is linear: half, 5 × 10⁻⁴ s (of whatever the 0.1 m change was).
8. L = ½ml²θ̇² + mgl cos θ. ∂L/∂θ̇ = ml²θ̇, d/dt → ml²θ̈; ∂L/∂θ = −mgl sin θ; Lagrange: ml²θ̈ + mgl sin θ = 0, i.e. θ̈ + (g/l) sin θ = 0.
