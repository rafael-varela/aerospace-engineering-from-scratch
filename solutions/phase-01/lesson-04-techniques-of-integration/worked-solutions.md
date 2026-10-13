# P01.04  Worked solutions

Seeded problems are worked for **s = 15**. Every number was produced by `solution.py`; if yours differs, run the code before trusting either. The lesson's wing: b = 11 m, s = 5.5 m, c_r = 1.6 m, τ = 0.6, σ = 12 kg/m², ℓ = 10.767 m.

## 1. Your own wing (s = 15)

b = 12 m, s = 6 m, c_r = 1.5 m, τ = 0.5, σ = 11.5 kg/m².

(a) Mass = σ c_r s (1 + τ) = 11.5 × 1.5 × 6 × 1.5 = **155.25 kg**.
(b) ȳ = s(1 + 2τ)/(3(1 + τ)) = 6 × 2/(3 × 1.5) = **2.667 m** (a rectangle would give 3 m, a triangle 2 m).
(c) I_xx = 2σ c_r s³[1/3 − (1 − τ)/4] = 2 × 11.5 × 1.5 × 216 × (0.33333 − 0.125) = 7452 × 0.20833 = **1552.5 kg m²**. `roll_inertia_numeric(lambda y: chord_linear(y, 1.5, 0.5, 6.0), 11.5, 6.0, 2, "simpson")` returns 1552.5000 to all digits because y² c(y) is a cubic and Simpson's error term carries f⁽⁴⁾ = 0.
(d) A rectangle of the same mass has chord c̄ = 155.25/(2 × 11.5 × 6) = 1.125 m and I = 2σ c̄ s³/3 = **1863.0 kg m²**, 20 % more: the taper moves mass inboard. Radius of gyration √(1552.5/155.25) = **3.162 m** = 0.527 s.

## 2. The four techniques

(a) Integration by parts, round one: u = y², dv = e^(−y/ℓ)dy, v = −ℓe^(−y/ℓ):
∫₀ˢ y² e^(−y/ℓ)dy = [−ℓ y² e^(−y/ℓ)]₀ˢ + 2ℓ ∫₀ˢ y e^(−y/ℓ)dy = −ℓ s² e^(−s/ℓ) + 2ℓ ∫₀ˢ y e^(−y/ℓ)dy.
Round two: u = y, same dv: ∫₀ˢ y e^(−y/ℓ)dy = [−ℓ y e^(−y/ℓ)]₀ˢ + ℓ ∫₀ˢ e^(−y/ℓ)dy = −ℓ s e^(−s/ℓ) + ℓ²(1 − e^(−s/ℓ)).
Combine: −ℓs²e^(−s/ℓ) − 2ℓ²s e^(−s/ℓ) + 2ℓ³ − 2ℓ³e^(−s/ℓ) = 2ℓ³ − e^(−s/ℓ)(s²ℓ + 2sℓ² + 2ℓ³). Multiply by 2σc_r: equation (3). The boundary terms at y = 0 vanish because of the factors y² and y; the only survivor of the lower limit is the 2ℓ³ from the last plain exponential.

(b) y = s sin θ, dy = s cos θ dθ, 0 ≤ θ ≤ π/2 covers 0 ≤ y ≤ s once; there cos θ ≥ 0 so √(1 − sin²θ) = +cos θ (the positive branch is the one the chord needs). The integral is s³ ∫₀^(π/2) sin²θ cos²θ dθ. Double angle: sin θ cos θ = ½ sin 2θ, so the integrand is ¼ sin²2θ = ⅛(1 − cos 4θ). ∫₀^(π/2) ⅛(1 − cos 4θ)dθ = ⅛[θ − sin 4θ/4]₀^(π/2) = π/16. So I_xx = 2σc_r s³ π/16 = πσc_r s³/8: equation (2).

(c) T − kv² = (√T − √k v)(√T + √k v). Set 1/(T − kv²) = A/(√T − √k v) + B/(√T + √k v); clearing denominators, 1 = A(√T + √k v) + B(√T − √k v); at v = √T/√k: A = 1/(2√T); at v = −√T/√k: B = 1/(2√T). Then
t = m ∫₀ᵛ dv′/(T − kv′²) = (m/(2√T)) [−(1/√k) ln(√T − √k v′) + (1/√k) ln(√T + √k v′)]₀ᵛ = (m/(2√(Tk))) ln[(√T + √k v)/(√T − √k v)],
and dividing numerator and denominator by √k gives (v_t + v)/(v_t − v). Expansion with x = v/v_t: ln((1 + x)/(1 − x)) = 2x + 2x³/3 + O(x⁵), so t ≈ (m/(2√(Tk)))(2v/v_t + 2v³/(3v_t³)) = mv/T + mkv³/(3T²) + …, using √(Tk) v_t = T and √(Tk) v_t³ = T²/k. The first term is the drag-free time; the correction is the extra time drag costs and grows as v³. Numbers (m = 1000 kg, T = 2500 N, k = 0.3675 kg/m): at 10 m/s the exact 4.0198 s against the two-term 4.0196; at 20 m/s 8.1626 against 8.1568; at 30 m/s 12.5756 against 12.5292 (x = 0.36, the x⁵ term is now 0.4 %).

(d) With h = (b − a)/n and Δf′ = f′(b) − f′(a): M_n = I − (h²/24)Δf′ + O(h⁴) and T_n = I + (h²/12)Δf′ + O(h⁴) (the lesson's §2.6, both from the Euler–Maclaurin expansion). Then (2M_n + T_n)/3 = I + (h²/3)(−2/24 + 1/12)Δf′ + O(h⁴) = I + O(h⁴): the h² terms cancel identically. The O(h⁴) coefficient of each rule is proportional to f‴(b) − f‴(a); the combination's is −(h⁴/180)[f‴(b) − f‴(a)], which is the Simpson error law, and it is zero whenever f‴ is constant, that is for every cubic. (2M_n + T_n)/3 *is* Simpson's rule with 2n panels: the midpoints become the odd nodes with weight 4h/6 and the trapezoid nodes the even ones.

## 3. Laplace transforms as improper integrals

(a) f = 1: ∫₀^∞ e^(−st)dt = lim [−e^(−st)/s]₀ᴿ = 1/s, for s > 0 (for s ≤ 0 the exponential does not decay).
f = t: parts with u = t, dv = e^(−st)dt: [−te^(−st)/s]₀ᴿ + (1/s)∫e^(−st)dt → 0 + 1/s² = 1/s², s > 0 (R e^(−sR) → 0 for s > 0).
f = e^(−at): ∫e^(−(s+a)t)dt = 1/(s + a), for s > −a.
f = sin ωt: parts twice with u = sin ωt then u = cos ωt: I = [−e^(−st) sin ωt/s]₀^∞ + (ω/s)∫e^(−st) cos ωt dt = (ω/s)([−e^(−st) cos ωt/s]₀^∞ − (ω/s)I) = (ω/s)(1/s − (ω/s)I). Solve: I(1 + ω²/s²) = ω/s², so I = ω/(s² + ω²), s > 0.

(b) With `improper_integral(lambda t: np.exp(-2*t)*f(t), 0.0, np.array([2, 5, 10, 20.]), n=4000)`:

| f | R = 2 | 5 | 10 | 20 | exact |
|---|---|---|---|---|---|
| 1 | 0.49084 | 0.499977 | 0.5000000 | 0.5 | 1/2 |
| t | 0.22711 | 0.249875 | 0.2500000 | 0.25 | 1/4 |
| e^(−t) | 0.33251 | 0.333333 | 0.3333333 | 0.3333333 | 1/3 |
| sin 3t | 0.22750 | 0.230773 | 0.2307692 | 0.2307692 | 3/13 = 0.2307692 |

All four reach the closed form by R = 10 (e^(−20) = 2 × 10⁻⁹). The sin case is not monotone (its tail oscillates), which is why the test in the lesson insists on monotonicity only for a positive integrand.

(c) f = eᵗ with s = 2: the integrand is e^(−t), and the sequence 0.86466, 0.99326, 0.99995, 1.00000 converges to 1 = 1/(s − 1). With s = 0.5 the integrand is e^(+0.5t) and the sequence is 3.44, 22.4, 295, 44 051: each step multiplies by about e^(0.5 ΔR). The transform of eᵗ exists only for s > 1; the general rule is s > (growth rate of f).

(d) A system whose response grows like e^(at) has a transform only for Re s > a; if you evaluate F(s) outside that region (for instance on the imaginary axis to get a frequency response) you are quoting a number for a divergent integral, and the frequency response you compute does not exist physically because the system is unstable.

## 4. Rules, bounds and the tip

(a) Errors |I_n − 1254.4379| on the raw elliptic integrand:

| n | midpoint | trapezoid | Simpson |
|---|---|---|---|
| 4 | 57.6 | 218 | 103 |
| 8 | 22.4 | 80.1 | 34.2 |
| 16 | 8.26 | 28.8 | 11.8 |
| 32 | 2.98 | 10.3 | 4.10 |
| 64 | 1.06 | 3.65 | 1.44 |
| 128 | 0.378 | 1.29 | 0.508 |
| 256 | 0.134 | 0.458 | 0.179 |

Log-log slopes on n ≥ 16: midpoint **1.49**, trapezoid **1.49**, Simpson **1.51**. Near the tip c ∝ √(s − y), so f is not differentiable there; every rule's error on the last panel is of order h × h^(1/2) = h^(3/2), independent of the weights, and that term dominates the h² (or h⁴) error from the smooth interior. Simpson's extra order needs four bounded derivatives, which the square root does not have.

(b) On the θ-form, 2σc_r s³ sin²θ cos²θ on [0, π/2], with n = 4: midpoint error 2 × 10⁻¹³, trapezoid 2 × 10⁻¹³, Simpson 2 × 10⁻¹³; all three are exact at n = 4. The distinction appears at n = 2: midpoint and trapezoid are still exact (error 10⁻¹³) and Simpson is off by 418 kg m². The integrand is (1 − cos 4θ)/8; on [0, π/2] the term cos 4θ is the *first* harmonic of the interval's period P = π/2. An equal-weight rule with n equally spaced nodes over a full period integrates the k-th harmonic exactly unless k is a non-zero multiple of n (only then do all nodes sample it at the same phase), so for k = 1 both the trapezoid and the midpoint rule are exact for every n ≥ 2: at n = 2 the trapezoid nodes give cos 4θ = 1, −1, 1 with weights ½, 1, ½, sum zero; the midpoint nodes give 0, 0. Simpson's 1-4-2-4-1 weights are not equal, so its n = 2 sum of cos 4θ is (h/3)(1 − 4 + 1) ≠ 0. The trapezoid rule is the right tool for periodic integrands, and its accuracy there is spectral, not second order.

(c) Exponential wing, exact 1459.139. max|f″| = f″(0) = 2 × 2σc_r = 76.8 (f = 2σc_r y² e^(−y/ℓ); numerically 76.74 on a fine grid), max|f⁽⁴⁾| = f⁽⁴⁾(0) = 2σc_r × 12/ℓ² = 3.975.

| rule | n | error | bound | bound/error |
|---|---|---|---|---|
| trapezoid | 8 | 7.428 | 16.64 | 2.24 |
| | 16 | 1.858 | 4.159 | 2.24 |
| | 32 | 0.4645 | 1.040 | 2.24 |
| | 64 | 0.1161 | 0.260 | 2.24 |
| | 128 | 0.02903 | 0.0650 | 2.24 |
| Simpson | 8 | 0.01805 | 0.02713 | 1.50 |
| | 16 | 0.001129 | 0.001696 | 1.50 |
| | 32 | 7.058 × 10⁻⁵ | 1.060 × 10⁻⁴ | 1.50 |

The ratio is constant because both error and bound scale as h² (or h⁴); the bound uses max|f″| over [0, s] while the true error is (h²/12)[f′(s) − f′(0)] = (h²/12) ∫f″ dy, the *average* of f″ times (b − a). The ratio 2.24 is max|f″| × s/∫f″ dy, a property of the integrand alone.

(d) Four significant figures means an error below 10⁻⁴ × 1459 = 0.146 kg m². Trapezoid: s³ × 76.8/(12n²) ≤ 0.146 → n ≥ 85.5 → **n = 86** (measured error 0.064). Simpson: s⁵ × 3.975/(180n⁴) ≤ 0.146 → n ≥ 5.5 → **n = 6** (measured 0.057). Simpson needs fourteen times fewer evaluations for the same guarantee on this smooth integrand; on the raw elliptic integrand the same guarantee would need n ≈ 1100 for Simpson, and four θ-panels.

## 5. Root bending moment of an elliptic load

(a) ∫₀ˢ √(1 − y²/s²) dy with y = s sin θ is s ∫₀^(π/2) cos²θ dθ = s π/4 (quarter of the ellipse's area). So L₀ πs/4 = W/2 and **L₀ = 2W/(πs) = 2 × 10 000/(π × 5.5) = 1157.5 N/m**.

(b) M₀ = ∫₀ˢ y L′ dy = L₀ ∫₀ˢ y √(1 − y²/s²) dy = L₀ s² ∫₀^(π/2) sin θ cos²θ dθ = L₀ s² [−cos³θ/3]₀^(π/2) = L₀ s²/3. With L₀ = 2W/(πs): M₀ = 2Ws/(3π) = **11 671 N m**. (The centroid of a quarter-ellipse's area is at 4s/(3π) = 2.334 m from the root, and M₀ = (W/2) × 2.334.)

(c) Load ∝ tapered chord: M₀ = (W/2) ȳ = 5000 × 2.521 = **12 604 N m**; uniform load: 5000 × 2.75 = **13 750 N m**. The elliptic distribution is kindest to the spar (8 % below the tapered, 15 % below the uniform) because it puts its load furthest inboard; it is also the minimum-induced-drag distribution (Phase 05), which is why sailplane designers chase it from both sides.

(d) Shear V(2) = ∫₂ˢ L′ dy and moment M(2) = ∫₂ˢ (y − 2)L′ dy, Simpson:

| n | V(2) N | M(2) N m |
|---|---|---|
| 100 | 2736.73 | 3959.18 |
| 200 | 2736.97 | 3960.02 |
| 400 | 2737.05 | 3960.32 |
| 1000 | 2737.09 | 3960.44 |
| exact | 2737.10 | 3960.48 |

The exact values come from the same substitution with θ₀ = arcsin(2/s). The errors fall by 2.8 per doubling (order 1.5, the tip again), so from the n = 200 → 400 change of 0.08 N the remaining error at n = 400 is about 0.08/(2^1.5 − 1) = 0.05 N: **V(2) = 2737.1 ± 0.1 N, M(2) = 3960.4 ± 0.2 N m**. The lesson's substitution would make n = 8 exact.

## 6. Where the model and the rules stop

(a) **Constant σ fails.** Two 60 kg point masses at 3.5 m add 2 × 60 × 3.5² = 1470 kg m², equal to the whole structural wing (1490.7): I_xx doubles with full tanks. Observable: the roll time constant with full and empty tanks differs by a factor of about 1.6 (the fuselage's share dilutes the ratio); a pilot feels it as a sluggish roll on departure.
(b) **Idealised planform.** Replacing the tapered chord by c_r over the inner 1.5 m adds only 2σ ∫₀^1.5 y²(c_r − c(y)) dy = 3.5 kg m², 0.2 %, because y² is small there; the fuselage (say 400 kg, radius 0.6 m, a thin shell) adds m r² ≈ 144 kg m² or, as a solid cylinder, m r²/2 = 72 kg m², 5–10 %. Observable: the measured I_xx exceeds the wing estimate by the fuselage's contribution, and the pendulum test (lesson part 6) can be done with and without the wings on some kit aircraft.
(c) **Tabulated data with a hidden singularity.** Twelve chords at 0.5 m spacing cover 5.5 m (n = 11 panels, odd, so Simpson is out); the trapezoid rule on the last panel treats √(h − u) as linear and gets (h/2)√h = 0.177 h^(3/2) instead of (2/3)h^(3/2) = 0.236: 25 % low on that panel. For the elliptic wing the total I_xx from the trapezoid rule on this table is 4.0 % low. Observable: refine the table near the tip (ask for chords at 5.25, 5.4, 5.45 m) and watch the answer rise; or fit c² (which is a quadratic for an ellipse) and integrate analytically.
(d) **The take-off model.** m dv/dt = T(v) − kv² − μ(W − L(v)) with T = 2500 − 26.7 v and L = ½ρv² S C_L (C_L ≈ 0.5 rolling, S = 15 m²). At v = 0 the acceleration is 1.715 m/s² (friction eats 0.785 of 2.5); at 30 m/s it is 0.915. Simpson on dv/a(v) gives t(30) = **23.6 s** against 12.58 s with drag only: the roll nearly doubles. The denominator is still a quadratic in v, so partial fractions still work with new roots. Observable: the take-off distance, which pilots measure against runway markers; on grass it is published as 15–30 % longer at the same weight, which says the real friction coefficient is nearer 0.04 than 0.08.

## 7. A heavier-at-the-root mass law

(a) Linear taper, u = y/s: I_xx = 2σ′ ∫₀ˢ y² c_r²[1 − (1 − τ)y/s]² dy = 2σ′c_r² s³ ∫₀¹ u²[1 − 2(1 − τ)u + (1 − τ)²u²] du = **2σ′c_r² s³ [1/3 − (1 − τ)/2 + (1 − τ)²/5]**. Mass: 2σ′c_r² s ∫₀¹ [1 − (1 − τ)u]² du = 2σ′c_r² s [1 − (1 − τ) + (1 − τ)²/3].
Elliptic: c² = c_r²(1 − y²/s²), so I_xx = 2σ′c_r² ∫₀ˢ (y² − y⁴/s²) dy = 2σ′c_r² s³(1/3 − 1/5) = **4σ′c_r² s³/15**; mass 2σ′c_r² s (1 − 1/3) = 4σ′c_r² s/3. No substitution: squaring removed the root.

(b) σ′ = 12/1.6 = 7.5 kg/m³, 2σ′c_r² = 38.4 (the same prefactor as before, by construction). Linear: bracket 0.33333 − 0.2 + 0.032 = 0.16533, **I_xx = 1056.3 kg m²**, mass 38.4 × 5.5 × 0.65333 = **138.0 kg**, k = √(I/m) = **2.767 m**. Elliptic: **I_xx = 851.8 kg m²**, mass **140.8 kg**, **k = 2.460 m**. Against λ ∝ c: inertias down 29 % and 32 % for masses down only 18 % and 15 %; the radii of gyration fall from 2.970 m and 2.750 m. The c² law puts the mass where the chord is big, which is inboard, and I_xx weights y². The elliptic wing is now 19 % below the tapered one instead of 16 %.

(c) The integrand y² c² is a quartic; Simpson's error law needs f⁽⁴⁾ = 0 for exactness, and here f⁽⁴⁾ = 2σ′ × 24 c_r²(1 − τ)²/s² = 0.9748, a constant. Simpson with n = 2 gives 1064.80, **8.518 kg m² (0.8 %) high**, and that is *exactly* the bound s⁵ f⁽⁴⁾/(180 × 2⁴) = 8.5184, because for constant f⁽⁴⁾ the error formula is an equality. The error falls as 16/n⁴: n = 4 gives 0.5324, n = 8 gives 0.0333, n = 32 gives 1.3 × 10⁻⁴. For 10⁻⁶ absolute: 136.29/n⁴ ≤ 10⁻⁶ → n ≥ 108.05 → **n = 110** (n = 108 gives 1.002 × 10⁻⁶). Since the error is known exactly, one Richardson step (16 S_{2n} − S_n)/15 makes n = 4 exact; that is Boole's rule.

(d) (i) Foam-core model wing: **λ ∝ c²**, because a solid section of constant thickness ratio has area ∝ c² and uniform density. (ii) Metal single-spar wing: **neither**; the spar cap area is sized by root bending moment and falls faster than c² outboard, while the skin goes as c; a λ = a c + b c² fit to a weight breakdown, or heavier than c² at the root. (iii) Sailplane with water ballast: **λ ∝ c² for the structure plus the ballast as a line mass** over the tank span; the ballast (often 100 kg or more) dominates I_xx exactly as in Problem 6(a), which is why ballasted sailplanes roll slowly.
