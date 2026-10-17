# P03.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. One step by hand (s = 15: h = 0.25)

**Derivation.** f(t, y) = −y + t; y(0.25) = 0.25 − 1 + 2e^{−0.25} = **0.807601566**.

(a) Euler: y₁ = 1 + 0.25 × (−1) = **0.75**, error −5.76 × 10⁻². Heun: k₁ = −1; predictor 1 + 0.25(−1) = 0.75; k₂ = −0.75 + 0.25 = −0.5; y₁ = 1 + 0.125(−1 − 0.5) = **0.8125**, error +4.90 × 10⁻³. RK4: k₁ = −1; k₂ = f(0.125, 0.875) = −0.75; k₃ = f(0.125, 1 + 0.125(−0.75)) = f(0.125, 0.90625) = −0.78125; k₄ = f(0.25, 1 + 0.25(−0.78125)) = f(0.25, 0.8046875) = −0.5546875; y₁ = 1 + (0.25/6)(−1 − 1.5 − 1.5625 − 0.5546875) = 1 − 0.1923828125 = **0.8076171875**, error +1.56 × 10⁻⁵. The three errors scale as h², h³, h⁵ with constants ½, 1/12, 1/120 (y″ = y′′ of the exact solution at 0 is 2, so Euler's h²y″/2 = 0.0625 ✓ to 8 %).

(b) Errors at t = 2:

| h | Euler | Heun | RK4 |
|---|---|---|---|
| 0.25 | 7.04 × 10⁻² | 6.89 × 10⁻³ | 2.17 × 10⁻⁵ |
| 0.125 | 3.45 × 10⁻² | 1.55 × 10⁻³ | 1.22 × 10⁻⁶ |
| 0.0625 | 1.71 × 10⁻² | 3.70 × 10⁻⁴ | 7.25 × 10⁻⁸ |

Ratios 2.04, 2.02 (Euler); 4.43, 4.20 (Heun); 17.8, 16.9 (RK4). Orders **1.02, 2.11, 4.11**.

(c) The error is Chᵖ + Dh^{p+1} + …, so the ratio is 2ᵖ(1 + (D/C)h)/(1 + (D/C)h/2), above 2ᵖ when D/C > 0 and tending to 2ᵖ as h → 0; the second ratio is closer than the first for all three methods ✓.

**Check.** Euler's global error ≈ (h/2)∫|y″|e^{−(2−t)}dt is of order h × 0.3: 0.07 at h = 0.25 ✓.

**Common mistakes.** Evaluating k₂ at the old y (giving Euler twice); forgetting the t argument advances with cᵢh; reading an order off a single ratio.

## 2. Local and global error, stability function, RK4 as Simpson

(a) From the exact solution, y(tₙ₊₁) = y(tₙ) + hf(tₙ, y(tₙ)) + h²y″(ξ)/2, and the Euler step is the first two terms: **local error h²y″(ξ)/2**. Subtracting the step from the exact recursion: eₙ₊₁ ≤ eₙ + h|f(tₙ, yₙ) − f(tₙ, y(tₙ))| + h²M/2 ≤ (1 + hL)eₙ + h²M/2 with M = max|y″|. With e₀ = 0, eₙ ≤ (h²M/2)Σ_{k<n}(1 + hL)ᵏ = (h²M/2)[(1 + hL)ⁿ − 1]/(hL) = (hM/2L)[(1 + hL)ⁿ − 1] ≤ **(hM/2L)(e^{LT} − 1)** since (1 + hL)ⁿ ≤ e^{nhL} = e^{LT}. ∎

(b) Heun on f = λy: k₁ = λy, k₂ = λ(y + hλy) = λy(1 + z), y₁ = y + (h/2)λy(2 + z) = **y(1 + z + z²/2)**. RK4: k₁ = λy; k₂ = λy(1 + z/2); k₃ = λy(1 + z/2 + z²/4); k₄ = λy(1 + z + z²/2 + z³/4); y₁ = y + (z/6)[1 + 2(1 + z/2) + 2(1 + z/2 + z²/4) + (1 + z + z²/2 + z³/4)]y = y + (z/6)[6 + 3z + z² + z³/4]y = **y(1 + z + z²/2 + z³/6 + z⁴/24)**. Heun's limit: on z = −x, R = 1 − x + x²/2 ≥ ½ > −1 always, and R ≤ 1 iff x²/2 ≤ x iff x ≤ 2: **limit 2**. RK4: R(−x) = 1 at x = 0 and crosses −1 at **x = 2.7853** (`stability_limit_real(4)`; between, |R| < 1). Euler on y′ = iωy: |R| = |1 + ihω| = √(1 + h²ω²) > 1 for every h > 0: unstable. ∎

(c) With f = g(t) the k's do not depend on y: k₁ = g(t), k₂ = k₃ = g(t + h/2), k₄ = g(t + h), and y₁ − y = (h/6)[g(t) + 4g(t + h/2) + g(t + h)], **Simpson's rule** on [t, t + h], exact for cubic g (`P03.04` 2.3). ∎

(d) The integration takes n = T/h steps, each adding a local error ≤ Ch^{p+1} that is then propagated with a factor at most e^{L(T−tₙ)}; the sum is ≤ e^{LT}C(T/h)h^{p+1} = O(hᵖ). The factor is unavoidable because the ODE itself may amplify perturbations (y′ = Ly does so exactly), but for y′ = −ky the propagation factor is e^{−k(T−tₙ)} < 1: early errors are forgotten and the global error is about C T hᵖ or less.

**Check.** RK4's local error on y′ = y at h = 0.5 is 2.84 × 10⁻⁴ against h⁵/120 = 2.60 × 10⁻⁴ (the h⁶/720 term adds 9 %) ✓.

**Common mistakes.** Claiming Euler's global error is h²; mixing z = hλ with h; forgetting that stability needs |R| ≤ 1, not R ≤ 1 (the lower crossing R = −1 is the one that bites RK4).

## 3. Apogee of a sounding rocket (s = 15: v₀ = 275 m/s)

**Assumptions.** m = 20 kg, ρ = 1.225 kg/m³, C_DA = 0.02 m², g = 9.80665 m/s², constant ρ (the real density falls 13 % over 1.4 km). SI units.

(a) v_t = √(2 × 20 × 9.80665/(1.225 × 0.02)) = **126.53 m/s**; t_a = (126.53/9.80665) arctan(275/126.53) = 12.903 × 1.1396 = **14.704 s**; z_a = (126.53²/(2 × 9.80665)) ln(1 + (275/126.53)²) = 816.3 × 1.7446 = **1424.13 m**. No drag: v₀²/2g = 3855.8 m, 2.7 times higher, and t = v₀/g = 28.0 s: drag at 2.2 v_t dominates the early coast.

(b) Altitude error at t = 14 s:

| h (s) | Euler | Heun | RK4 |
|---|---|---|---|
| 0.5 | 3.53 | 0.620 | 1.67 × 10⁻³ |
| 0.25 | 2.06 | 0.139 | 9.62 × 10⁻⁵ |
| 0.125 | 1.09 | 0.0331 | 5.78 × 10⁻⁶ |

Orders **0.84, 2.11, 4.09** (Euler ratios 1.7 and 1.9: still approaching 2).

(c) Apogee time error (s) / apogee height error (m):

| h (s) | Euler | Heun | RK4 |
|---|---|---|---|
| 0.5 | −0.451 / +1.41 | +0.0066 / −0.575 | +2.8 × 10⁻⁵ / −1.63 × 10⁻³ |
| 0.25 | −0.222 / +1.07 | +0.0016 / −0.128 | +3.9 × 10⁻⁶ / −8.4 × 10⁻⁵ |
| 0.125 | −0.110 / +0.702 | +0.0004 / −0.030 | +5.8 × 10⁻⁷ / −5.4 × 10⁻⁶ |

Heun's errors fall by 4 and RK4's by 16 per halving (orders 2.1 and 4.1 for z_a; 2.0 and 2.8–2.9 for t_a, where the linear interpolation of v adds its own O(h²) error that RK4 cannot beat). Euler's apogee time converges at order 1 (ratios 2.0) but its height error only at ratios 1.3 and 1.5: the Euler altitude is too high (the velocity error is positive: Euler under-estimates drag, which acts on the *old* velocity) and the parabola fit over a 0.5 s bracket adds an error of opposite sign that cancels part of it at the coarse steps. Verdict: Heun and RK4 converge at their orders in z_a; Euler is not yet asymptotic.

(d) 0.1 m: RK4 at h = 0.5 s is already at 1.6 mm (and 1.7 mm at t = 14 s); even h = 1 s gives 3.1 cm at 14 s. **RK4 with h = 1 s**: 15 steps, **60 evaluations** for a 15 s prediction; Heun would need h ≈ 0.2 s (150 evaluations) and Euler h ≈ 0.01 s (1500). The real limit is the model (constant ρ, constant C_D), not the integrator.

**Check.** v̇(0) = −9.81 − (1.225 × 0.02/40) × 275² = −9.81 − 46.3 = −56.1 m/s²; Euler's first step over-predicts v by about ½ v̈ h² with v̈ = 2(ρC_DA/2m)v|v̇| ≈ 2 × 6.1 × 10⁻⁴ × 275 × 56 = 18.9 m/s³: 2.4 m/s at h = 0.5, which after 14 s of coast is metres of altitude ✓ (3.5 m observed).

**Common mistakes.** Writing drag as −kv² (wrong sign after apogee); detecting apogee as the grid maximum of z (an error of ½|v̇|(h/2)² ≈ 0.3 m at h = 0.5 even for an exact integrator); fitting an order to Euler's three apogee values and reporting 0.5.

## 4. A different orbit (s = 15: a = 11 000 km, e = 0.25)

(a) T = 2π√(a³/μ) = **11 481.5 s** (3.19 h); perigee altitude a(1 − e) − R_E = 8250 − 6371 = **1879 km**; v_p = √(μ(1 + e)/(a(1 − e))) = **7771.4 m/s**, v_a = √(μ(1 − e)/(a(1 + e))) = **4662.8 m/s** (ratio 5/3 = (1 + e)/(1 − e) ✓).

(b) One-orbit position error (m):

| steps/orbit | Heun | Kutta3 | RK4 |
|---|---|---|---|
| 200 | 2.18 × 10⁵ | 5.04 × 10³ | 7.65 |
| 400 | 5.27 × 10⁴ | 6.30 × 10² | 0.421 |
| 800 | 1.30 × 10⁴ | 78.8 | 0.0245 |

Orders **2.04, 3.00, 4.14**.

(c) h = (T/800)(1/0.0245)^{1/4} = 2.53 T/800 → **317 steps per orbit**; a run gives **1.10 m** ✓.

(d) RK4 at 800 steps: 0.0245, 0.0525, 0.119 m after 1, 2, 4 orbits: growth factors 2.14 and 2.27 per doubling, i.e. **linear in time** (slightly faster while the error vector turns along-track). One day is 86 400/11 481.5 = 7.525 orbits, so the 800-step error after a day is about 7.5 × 0.0245 × 1.1 ≈ 0.20 m; (2) with p = 4 gives h = (T/800)(1/0.20)^{1/4} → **525 steps per orbit**; a run gives **2.08 m**, not 1 m: the one-day error at 525 steps grows faster than the pure h⁴ law assumed (order over this range is 4.14, and the linear-growth factor was taken from the first four orbits). With a factor-2 margin, 650 steps per orbit (**4900 steps, 19 600 evaluations** for the day) is the honest answer.

(e) After five orbits at 400 steps per orbit: RK4 **−4.2 × 10⁻⁹** (dissipative, −8.5 × 10⁻¹⁰ per orbit), Heun **+1.4 × 10⁻⁴** (anti-dissipative). The signs are the lesson's 2.7.

**Check.** RK4's 400-step error at e = 0.25 (0.42 m) is 5.3 times the lesson's e = 0.1 value (0.079 m) despite the longer period: the perigee speed ratio penalises a fixed step ✓ (lesson's extension 5).

**Common mistakes.** Using the one-orbit error for the day without multiplying by the orbit count; trusting (2) without a margin when the measured order is not exactly 4; reporting the energy drift without its sign.

## 5. Stiffness and the stability limit (s = 15: λ = 250)

(a) h_max = 2/250 = **0.00800 s** for Euler and Heun, 2.7853/250 = **0.01114 s** for RK4. Runs to t = 10 (exact y(10) = cos 10 = −0.8391 to 4 × 10⁻³): at 0.95 h_max, |y(10)| = 0.8412, 0.8410, 0.8411 (Euler, Heun, RK4); at 1.05 h_max, **7 × 10⁴⁸, 1.5 × 10⁵¹, 7.6 × 10⁷⁶**. The limit is exactly where the theory puts it.

(b) Steps to t = 10 at 0.95 h_max: Euler 1316 (1316 evaluations), Heun 1316 (2632), RK4 945 (3780). RK4's error there is **9.8 × 10⁻⁵**. Accuracy alone: errors 3.96 × 10⁻⁵ at h = 0.01 and 6.2 × 10⁻⁷ at 0.005 (ratio 64: the transient's h⁵ term is still in play), and (2) with p = 4 from the 0.005 run gives h ≈ 0.0056 for 10⁻⁶. So accuracy asks for h ≈ 0.0056 and stability for 0.0106: here the two happen to be close; for a smoother target (λ = 2500 or a tolerance of 10⁻⁴) stability would dominate by a factor 10–100. The problem is mildly stiff.

(c) `solve_ivp` at rtol = 10⁻⁶: **RK45 8138 evaluations**, **Radau 703**, **LSODA 495**, all with errors near 1.5 × 10⁻⁷. The explicit RK45 is hostage to the stability limit (its step is pinned near 3.3/λ for the whole run, 8000 evaluations for a solution that is cos t); the implicit Radau and the stiff-switching LSODA take steps of 0.1–0.5 once the transient is over. This is `P03.06`.

(d) Euler on y′ = iy: |R(ih)|² = 1 + h²; after n = 2π/h steps the amplitude grows by (1 + h²)^{π/h}: **1.368 at h = 0.1, 1.032 at h = 0.01** per period (≈ e^{πh}). RK4: |R(ih)| = 0.99999999306 at h = 0.1 and 1 − 7 × 10⁻¹⁵ at h = 0.01; per period 0.99999956 and 1 − 4 × 10⁻¹². Prediction: an Euler orbit spirals outward by 37 % in radius per revolution at 63 steps per orbit and by 3 % at 628, whatever the accuracy; an RK4 orbit decays by 4 × 10⁻⁷ per revolution at 63 steps (the lesson's −2.3 × 10⁻⁹ energy drift per ten orbits at 400 steps is this effect at h = 0.0157).

**Check.** 1.05 h_max for Euler gives |R| = |1 − 2.1| = 1.1 per step; 1250 steps: 1.1¹²⁵⁰ = 10⁵¹.7 — observed 7 × 10⁴⁸ over 1190 steps: 1.1¹¹⁹⁰ = 10⁴⁹.3 ✓.

**Common mistakes.** Using λ instead of −λ in `max_stable_step`; concluding from (b) that RK4 is "stable enough" (λ is small here by design); comparing nfev at different tolerances.

## 6. A method and a step for a one-day propagation

**Assumptions.** Two-body truth; 13.37 orbits per day; errors grow linearly in time once settled; machine rate 2 × 10⁶ evaluations/s; 2000 objects.

(a) RK4 at 800 steps per orbit: one-day error measured 0.143 m (13.4 × 4.7 mm × 2.3: the growth settles to about 6.6 mm per orbit). For 10 m, (2) gives h = (T/800)(10/0.143)^{1/4} = 2.89 T/800 → 277 steps per orbit, but a run gives 15.3 m (the order over this range is 4.17 and the day-long growth is not perfectly linear). With margin: **RK4 at 400 steps per orbit** gives **2.9 m** after a day; 330 gives 7.0 m. At 400: 5350 steps per object per day, **21 400 evaluations per object, 4.3 × 10⁷ in total, 21 s of machine time**. Recommendation: RK4, 400 steps per orbit (16 s step), with the Kepler check run as part of the tool's tests.

(b) Heun: one-day error at 800 steps ≈ 13.4 × 4.8 km = 64 km; 10 m needs h/80: 64 000 steps per orbit, **1.7 × 10⁶ evaluations per object**, 3.4 × 10⁹ in total, **29 minutes**. Euler: 4 × 10⁹ steps per orbit, 5 × 10¹⁰ per object-day, **10¹⁴ in total, 1.6 years**. Never.

(c) At e = 0.7 the perigee speed is (1 + e)/(1 − e) = 5.7 times the apogee speed; a step tuned to apogee is 5.7 times too long at perigee and the h⁴ law makes the error 10³ times larger there (a run at e = 0.7, a = 25 000 km, 245 steps per orbit: 22 km after one orbit). Tell the author: use an adaptive step controlled by a local error estimate (`P03.06`'s embedded RK45), or at minimum a step proportional to r^{3/2} (constant steps in eccentric anomaly), and verify against `kepler_state` for each eccentricity class in the catalogue.

(d) At rtol = 10⁻³ `solve_ivp` accepts a relative error of 10⁻³ in a state of magnitude 7.5 × 10⁶ m: the oracle itself is wrong by kilometres, 100 times the requirement, and its default RK45 can be *less* accurate than the fixed-step RK4 under test. The right oracle is the closed form (`kepler_state`, exact to rounding), and where no closed form exists, `solve_ivp` with `DOP853` at rtol = 10⁻¹³ (15 µm on this orbit), plus the self-consistency checks that need no oracle: the observed order and the energy drift.

**Check.** 21 400 evaluations per object-day for 2.9 m against 1.7 × 10⁶ for Heun: the factor 80 in cost is (10³)^{1/2}/(10³)^{1/4}-type arithmetic of the two orders ✓.

**Common mistakes.** Sizing the step from the one-orbit error without the growth in time; forgetting the margin for a measured order that is not exactly 4; validating an integrator against a less accurate one.

## 7. Order conditions of two-stage methods and the stability function

(a) For y′ = f(y): y″ = f′(y)y′ = f′f; y‴ = (f′f)′ = f″f·f + f′·f′f = f″f² + f′²f. Hence **y(t + h) = y + hf + h²f′f/2 + h³(f″f² + f′²f)/6 + O(h⁴)**. ∎

(b) k₁ = f. k₂ = f(y + c₂hf) = f + c₂hf′f + c₂²h²f″f²/2 + O(h³). Step: y + h(b₁ + b₂)f + h²b₂c₂f′f + h³b₂c₂²f″f²/2 + O(h⁴). ∎

(c) Matching h: **b₁ + b₂ = 1**; h²: **b₂c₂ = ½**. These are (1). At h³ the exact solution has f″f²/6 + f′²f/6; the method has b₂c₂²f″f²/2 and *no* f′²f term (it would need a stage evaluated at a point that itself depends on a previous slope, i.e. a third stage). Since f′²f is not identically zero (take f = y: f′²f = y), no two-stage explicit method is third order. ∎

(d) Induction: k₁ = λy = λyP₁(z) with P₁ = 1. If kⱼ = λyPⱼ(z) with deg Pⱼ ≤ j − 1 for j < i, then kᵢ = λ(y + hΣⱼaᵢⱼkⱼ) = λy(1 + Σⱼaᵢⱼ z Pⱼ(z)) = λyPᵢ(z) with deg Pᵢ ≤ i − 1. The step is y + hΣbᵢkᵢ = y(1 + zΣbᵢPᵢ(z)) = yR(z) with **deg R ≤ s**. The exact solution multiplies y by e^z; a method of order p has local error O(h^{p+1}) on every smooth problem, in particular on y′ = λy, so R(z) − e^z = O(z^{p+1}), i.e. R agrees with Σ_{k≤p}z^k/k! through degree p. When s = p the polynomial R has degree ≤ p and agrees with the truncated exponential through degree p: **R(z) = Σ_{k≤p}z^k/k!** exactly. For s = p = 1, 2, 3, 4 these are the functions of 2.5. ∎

**Check.** Problem 2(b) computed Heun's and RK4's R(z) by hand and found the truncated series ✓. For s > p (RK45's fifth-order method has s = 6) R has extra terms beyond z⁵ and the stability region differs from the truncated exponential's, which is why `P03.06` recomputes it.

**Common mistakes.** Treating f′ as ∂f/∂t; claiming order 3 from matching only the f″f² term; assuming R is the truncated exponential for every RK method regardless of stage count.

## Check-yourself answers

1. Euler local 2, global 1; Heun 3 and 2; RK4 5 and 4. There are T/h steps, so the global error is the local error divided by h.
2. b₁ + b₂ = 1 and b₂c₂ = ½; midpoint (c₂ = ½, b = (0, 1)) and Heun (c₂ = 1, b = (½, ½)).
3. R(z) = 1 + z + z²/2 + z³/6 + z⁴/24; real-axis limit |hλ| = 2.7853; for λ = −1000, h must be below 0.00279 or the solution grows without bound.
4. On an oscillator z = ±iωh and |1 + iωh| = √(1 + ω²h²) > 1 for every h: the amplitude grows every step, so a numerical orbit spirals outward and its energy rises.
5. h = h₀(1/0.0047)^{1/4} ≈ 3.8h₀: about 210 steps per orbit; it assumes the error is already in its h⁴ regime (check the observed order first).
6. Simpson's rule on [t, t + h]; exact for cubic g.
7. For the orbit's imaginary eigenvalues ±iω, |R(iθ)| for RK4 is slightly below 1 (1 − θ⁶/144 + …), so each step removes a little energy, at a constant rate per step: negative and linear. Euler has |R| = √(1 + θ²) > 1: energy added every step.
8. The stability limit 2.79/|λ_max| decides the step, not the accuracy; the remedy is an implicit method whose stability region contains the left half-plane (`P03.06`).
