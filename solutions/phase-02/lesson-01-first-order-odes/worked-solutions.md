# P02.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. A burn by hand ($s = 15$: ṁ = 6.5 kg/s)

**Assumptions.** Vertical, vacuum, constant g, ṁ and v_e. Units SI throughout.

**Derivation.** (a) t_b = (300 − 120)/6.5 = 27.69 s. Δv = 2400 ln(300/120) = 2400 × 0.91629 = 2199.1 m/s. Gravity loss g t_b = 9.81 × 27.69 = 271.7 m/s. v_b = 2199.1 − 271.7 = **1927.4 m/s**. (b) h_b = 2400 × 27.69 − (2400 × 120/6.5) ln 2.5 − ½ × 9.81 × 27.69² = 66 462 − 44 308 × 0.91629 − 3761 = 66 462 − 40 599 − 3761 = **22 101 m**. (c) Apogee = 22 101 + 1927.4²/19.62 = 22 101 + 189 348 = **211.4 km**.

**Result.** `vacuum_velocity(t_b, p)` = 1927.436 m/s, `vacuum_altitude` = 22 101.34 m; RK4 with 400 steps: 1927.436 m/s and 22 101.34 m, errors 9 × 10⁻¹⁰ m/s and 2 × 10⁻⁸ m (rounding-level; the truncation error at h = 0.07 s is below 10⁻⁹).

**Check.** Dimensions of each term of h_b are metres (m/s × s; (m/s × kg)/(kg/s); m/s² × s²). Limiting case: with g = 0 the apogee is infinite (no force to stop the rocket), and with v_e → 0 everything vanishes.

**Common mistakes.** Using m₀ instead of m(t_b) in the prefactor of the logarithm in (4); forgetting the gravity loss; using the lift-off T/W to judge the burnout acceleration.

## 2. The altitude integral and the drag coast

**Derivation.** (a) h(t) = ∫₀ᵗ v_e ln(m₀/(m₀ − ṁτ)) dτ − ½gt². Substitute u = m₀ − ṁτ, dτ = −du/ṁ, limits u: m₀ → m(t):

∫₀ᵗ ln(m₀/(m₀ − ṁτ)) dτ = (1/ṁ) ∫_{m(t)}^{m₀} ln(m₀/u) du.

Integrate by parts with f = ln(m₀/u), g′ = 1: ∫ ln(m₀/u) du = u ln(m₀/u) − ∫ u · (−1/u) du = u ln(m₀/u) + u. Evaluate: [u ln(m₀/u) + u]_{m(t)}^{m₀} = (0 + m₀) − (m ln(m₀/m) + m) = m₀ − m − m ln(m₀/m). Since m₀ − m = ṁt, the first integral is (v_e/ṁ)(ṁt − m ln(m₀/m)) = v_e t − (v_e m/ṁ) ln(m₀/m). Add −½gt². At t = 0, m = m₀, ln 1 = 0 and h(0) = 0 ✓.

(b) dh/dv = (dh/dt)/(dv/dt) = v/(−g − kv²). Δh = ∫_{v_b}^{0} −v dv/(g + kv²) = ∫₀^{v_b} v dv/(g + kv²). With w = g + kv², dw = 2kv dv: = (1/2k) ∫ dw/w = (1/2k) ln((g + kv_b²)/g) = **ln(1 + kv_b²/g)/(2k)**. The drag term was written −kv², which is −kv|v| only for v > 0: at the apex v changes sign, the drag term becomes +kv², and a different integral applies on the way down. The formula stops at the apex.

(c) With x = kv_b²/g: Δh = (1/2k)(x − x²/2 + x³/3 − …) = **v_b²/(2g) − k v_b⁴/(4g²) + k²v_b⁶/(6g³) − …**. The first term is the vacuum height; the second is the first-order drag loss, proportional to k and to v_b⁴ (drag ∝ v² acting over a distance ∝ v²). The series converges for x < 1 only. Numerically for v_b = 1000 m/s: k = 10⁻⁶ m⁻¹ gives x = 0.102, `quadratic_drag_coast` → 48 535 m, two-term expansion 48 371 m, difference 164 m, close to the third term 176 m ✓. For k = 10⁻⁵, x = 1.02: exact 35 139 m, two-term expansion 24 991 m, hopeless, because the series does not converge.

**Common mistakes.** Dropping the "+u" from the integration by parts; using the upward formula past the apex; trusting a Taylor expansion outside its radius of convergence.

## 3. Avionics box ($s = 15$: C = 2750 J/K)

**Assumptions.** Lumped thermal mass (uniform temperature), constant conductance and power. Units: J/K, W/K, W, K.

**Derivation.** (a) First order, linear with constant coefficients, autonomous (hence also separable). Standard form: dT/dt + (G/C) T = (P + GT_a)/C. (b) Integrating factor e^(Gt/C): d/dt[T e^(Gt/C)] = ((P + GT_a)/C) e^(Gt/C); integrate: T = T_∞ + (T₀ − T_∞) e^(−t/τ) with **τ = C/G = 1833 s = 30.6 min** and **T_∞ = T_a + P/G = 20 + 20 = 40 °C**. (c) 15 = 40 − 50 e^(−t/τ) → e^(−t/τ) = 0.5 → t = τ ln 2 = **1271 s = 21.2 min**. (d) RK4 with 200 steps to 1270.77 s gives T = 14.99999999998 °C (error 2 × 10⁻¹¹ K). `equilibria` on f(T) = (P − G(T − T_a))/C returns [(40.0, True)]: one stable equilibrium; the phase line's arrows all point at 40 °C, so the box warms monotonically toward T_∞ and never overshoots.

**Check.** After one time constant T = 40 − 50/e = 21.6 °C ✓ (63 % of the way). If P = 0 the equation is Newton's law of cooling toward T_a with the same τ.

**Common mistakes.** Forgetting that the steady state is *above* the cabin temperature because the box heats itself; quoting τ as "the time to reach steady state" (it is the time to cover 63 % of the gap).

## 4. Your rocket's drag ($s = 15$: C_D = 0.45)

**Result.** (a) RK4, 3000 steps to 150 s: burnout v_b = 1223.7 m/s, h_b = 11 831 m; apogee 51 049 m at t = 106.4 s. (b) Apogee with 300, 600, 1200, 3000 steps: 51 049.11, 51 049.12, 51 049.18, 51 049.18 m: converged to 0.1 m (the small non-monotone wobble is the grid resolution of `argmax`, not the integrator). (c) ρ(h_b) = 0.304 kg/m³, k = 4.84 × 10⁻⁵ m⁻¹: frozen-density coast gains 21 962 m in 55.9 s; vacuum coast gains 76 320 m in 124.7 s; the full model gains 39 218 m, between the two ✓. (d) `solve_ivp` at rtol 10⁻¹⁰: v_b = 1223.6850021 m/s, RK4 1223.6850021: difference 9 × 10⁻⁹ m/s.

**Check.** C_D = 0.45 lies between the lesson's 0.3 (apogee 75.1 km) and 0.5 (45.6 km) cases and so does the apogee, 51.0 km ✓.

**Common mistakes.** Reading the burnout state at index 20/Δt without checking `ts[i_b]`; using m₀ instead of m_dry in k.

## 5. The way down

**Result.** (a) Terminal velocities (100 kg, C_D = 0.5, A = 0.0707 m²): 45 km: 3004 m/s; 30 km: 1243; 20 km: 690; 10 km: 383; 0 km: 213 m/s. (b) Descent from the 45.56 km apogee (t = 99.8 s), impact at t = 214.1 s. At 30 km the rocket falls at 537 m/s, only 43 % of the local terminal speed (it is still accelerating from rest at apogee); at 18 km, 631 m/s against 614 (it has just crossed); at 10 km, 554 against 383 (44 % too fast); at the ground, **284 m/s against 213 m/s**, 33 % faster than terminal. The descent is *not* quasi-steady: v crosses the terminal curve once and then stays above it, because the air thickens faster than drag can slow the body. (c) The relaxation time toward the equilibrium is τ = 1/(2k|v*|) = 63 s at 30 km, 31 s at 18 km, 20 s at 10 km, 11 s at sea level, while the rocket crosses a scale height in 14–30 s. Where τ exceeds the crossing time the equilibrium moves before the velocity can reach it; only below about 3 km is τ comfortably shorter, and by then the lag accumulated above has not been shed. Impact at 284 m/s rather than 213 is the quantitative consequence. (d) The first assumption to fail is C_D = 0.5 nose-first: an unguided empty rocket body tumbles on the way down, and a tumbling cylinder has a C_D of order 1–1.5 referenced to a much larger area, so the real descent is slower and nearer quasi-steady. The measurable sign is a recovered body with impact damage on its side rather than its nose, or telemetry showing a descent rate well below the nose-first prediction.

**Common mistakes.** Assuming terminal velocity is reached "quickly"; computing v_t with the sea-level density at every altitude.

## 6. The global error of Euler's method

**Derivation.** (a) With f independent of y, Euler gives y_N = y₀ + h Σ_{n=0}^{N−1} f(t_n): the left Riemann sum L_N of ∫_a^b f dt. The trapezoid rule is T_N = L_N + ½h[f(b) − f(a)] and has error O(h²) (P01.03), so ∫f dt − L_N = ½h[f(b) − f(a)] + O(h²). The global Euler error is therefore ½h[f(b) − f(a)] to leading order. For the vacuum burn f = T/m − g runs from 56.2 to 155.2 m/s², so the error constant is ½(155.2 − 56.2) = 49.5 m/s per second of step; the lesson measures 47.8, 48.6, 49.1, 49.3 for h = 2, 1, 0.5, 0.25 s, converging to 49.5 ✓.

(b) Let e_n = |y(t_n) − y_n|. Taylor: y(t_{n+1}) = y(t_n) + h f(t_n, y(t_n)) + ½h²y″(ξ). Subtract the Euler step and use the Lipschitz bound |f(t, y) − f(t, y_n)| ≤ L|y − y_n|: e_{n+1} ≤ (1 + hL) e_n + ½h²M. Unroll with e₀ = 0: e_n ≤ ½h²M Σ_{j=0}^{n−1}(1 + hL)^j = ½h²M [(1 + hL)ⁿ − 1]/(hL) = (hM/2L)[(1 + hL)ⁿ − 1] ≤ (hM/2L)(e^{nhL} − 1) = **(hM/2L)(e^{L(t_n − t₀)} − 1)**, using 1 + x ≤ eˣ.

(c) k = 1.225 × 0.5 × 0.0707/(2 × 100) = 2.165 × 10⁻⁴ m⁻¹; v* = 212.9 m/s; L = 2k v* = 0.0922 s⁻¹; M = max|f′ f| = 2k v*(g + kv*²) = 1.808 m/s³. Over 60 s, e^{60L} − 1 = 251.3, so the bound is 2464 h m/s: 2464, 1232, 616, 308, 154 m/s for h = 1, 0.5, 0.25, 0.125, 0.0625 s. The measured errors (against `solve_ivp`, rtol 10⁻¹², reference v(60) = −211.196 m/s) are 0.307, 0.157, 0.079, 0.040, 0.020 m/s: first order as expected, and the bound is **about 8000 times pessimistic**. The reason is the sign of ∂f/∂v = −2k|v| < 0: the equation is *contracting*, errors made early decay instead of growing, while the bound assumes the worst case |∂f/∂y| = L with the growing sign. The bound is sharp only for expanding equations such as y′ = Ly.

**Common mistakes.** Taking L as the maximum of |f| instead of |∂f/∂y|; forgetting that the bound is for the worst case and concluding Euler is useless.

## Check-yourself answers

1. If f and ∂f/∂y are continuous near (t₀, y₀), exactly one solution passes through (t₀, y₀) on some interval. dy/dt = 3y^{2/3}, y(0) = 0 has y = 0 and y = t³: ∂f/∂y is unbounded at y = 0.
2. dm = −ṁ dt → m = m₀ − ṁt, valid for 0 ≤ t ≤ t_b = (m₀ − m_dry)/ṁ; after that the mass is constant and a different equation (no thrust) applies.
3. Because dv = v_e dm/m integrates to a logarithm: each kilogram of propellant gives more Δv when the rocket is lighter, so the gain depends on the ratio m₀/m_f, and doubling the propellant does not double Δv.
4. The term g t_b in v_b: the velocity gravity removes during the burn. It depends only on the burn duration, so a faster burn (higher ṁ for the same propellant) loses less to gravity, at the cost of higher drag in dense air.
5. A function μ(t) such that μ(y′ + p y) = (μ y)′, namely μ = e^{∫p dt}. For dv/dt + (c/m)v = −g, μ = e^{ct/m}.
6. At v* = −√(g/k), the terminal velocity. Stable, because f′(v*) = −2k|v*| < 0: a perturbation decays with time constant 1/(2k|v*|).
7. Local error: the error of one step from exact data, ½h²y″ for Euler. Global error: the accumulated error at a fixed final time after (b − a)/h steps, O(h) for Euler. RK4: local O(h⁵), global O(h⁴).
8. Drag only removes energy, so the vacuum coast is an upper bound; freezing the density at the burnout value overstates drag everywhere above burnout (the real density falls by e every 8.5 km), so it is a lower bound. The frozen-density assumption is the pessimistic one.
