# P03.04  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `p = EntryParams()`, h = (40 + s/2) × 10³ m. 2. `deceleration_analytic(h, p)` and `deceleration_from_velocity(h, p, dh)` for dh = 100, 50, 25; `observed_order`. 3. V‴ by the five-point third-difference stencil with h = 20 m (or symbolically: V = V_e exp(−Bu), u = e^{−h/H}; V′ = (B/H)uV, and differentiate twice more). Predicted error = V sin γ · dh² V‴/6. 4. h* from `optimal_step_central(V, |V‴|)`; sweep dh downwards and watch the error.

## 2.
(a) f(x ± h) = f ± hf′ + h²f″/2 ± h³f‴/6 + h⁴f⁗/24 ± …; forward: subtract f, divide by h; central: subtract the two expansions, divide by 2h; second: add, subtract 2f, divide by h². (b) The linear interpolation error on [0, h] is ½f″(ξ)t(t − h); ∫₀ʰ t(t − h)dt = −h³/6; multiply by ½; sum n panels with nh = b − a, replacing the ξ's by one by the intermediate value theorem. (c) Centre at 0: Simpson on [−h, h] is (h/3)(f(−h) + 4f(0) + f(h)); for f = x and x³ both sides are 0; for 1 and x² both sides agree; so exact to degree 3. For x⁴: exact ∫ = 2h⁵/5, Simpson = 2h⁵/3, error −4h⁵/15 = −(h⁵/90)f⁗ with f⁗ = 24; per pair of panels (width 2h) that is −(2h)h⁴f⁗/180; sum over n/2 pairs. (d) dV/dh = (dV/dt)/(dh/dt) = ρV/(2β sin γ); with ρ = ρ₀e^{−h/H}: d ln V = (ρ₀/(2β sin γ))e^{−h/H}dh; integrate from ∞ where V = V_e. a = −dV/dt = V sin γ dV/dh. Set d/dh[ρV²] = 0 using d ln V/dh = Bu/H and d ln ρ/dh = −1/H with u = e^{−h/H}: −1 + 2Bu = 0. Same for ½ ln ρ + 3 ln V: −½ + 3Bu = 0. a_max: substitute u = 1/2B. (e) Q = ∫q̇ dt = ∫q̇ dh/(V sin γ); q̇/V = k√(ρ₀/Rₙ)u^{1/2}V_e²e^{−2Bu}; dh = −H du/u; the integral ∫u^{−1/2}e^{−2Bu}du = √(π/2B) erf(√(2Bu)).

## 3.
(a) Triangle ½ × 2000 × 0.2 plus 2000τ(1 − e^{−5.8/τ}). (b) `t = numpy.linspace(0, 6, n + 1)`, F by `numpy.where`; `trapezoid`/`simpson` accept a callable, so pass F; n = 60 and n = 40. (c) n = 120. Predicted factors 4 (trapezoid) and 16 (Simpson). For the trapezoid rule compute the Euler–Maclaurin end term −h²[F′(6) − F′(0⁺)]/12 and the kink term h²[F′(0.2⁺) − F′(0.2⁻)]/12 (a kink on a node is an end of two sub-intervals); their sum changes sign with h because the kink term is h²-exact while the tail term is tiny, so the error is not monotone in h.

## 4.
(a) `p = EntryParams(beta=200 + 10 s)`; `entry_B`, `peak_heating_altitude`, `peak_deceleration`, `heat_load_analytic(120e3, p)`. Scaling: Q(β)/Q(300) against √(300/β); h_q(β) − h_q(300) against H ln(300/β). (b) Loop over methods and n; `observed_order` on n = 32, 64, 128 with h = 120 km/n. (c) Increase even n from 2 until the relative error drops below 10⁻⁶.

## 5.
(a) Reference by `scipy.integrate.quad`; loop. (b) `romberg(f, 0, 2, 6)`; the diagonal entry R[i, i] is a Newton–Cotes-like rule of degree 2i on 2ⁱ panels; its error constant involves f^(2i+2), which for e^{−x}cos 3x grows like 3^{2i+2}, so for coarse panels the higher rule can lose. (c) Exact value 1/(a + 1). Orders from `observed_order`; the trapezoid and Simpson orders tend to a + 1; Gauss converges algebraically with order about 2a + 2 because the nearest singularity is on the interval, not off it. `quad` subdivides adaptively and extrapolates (QUADPACK's qags), which recovers full accuracy.

## 6.
(a) Newton–Cotes rules (trapezoid, Simpson, Romberg over them) run on equally spaced samples; Gauss and adaptive methods need a callable. (b) Read Simpson's 10⁻⁶ crossing from the n table (120 km / n). (c) Simpson on all samples S_h and on every other sample S_{2h}; estimate |S_h − S_{2h}|/15; evaluate with `heat_load_numeric(120e3, p, "simpson", 120)` and `(…, 60)`. (d) Separate numerical error (controllable, cheap) from model error (not), and say which one the requirement should be about.

## 7.
(a) The rule integrates the interpolant exactly, and the interpolant of a degree ≤ n − 1 polynomial is itself. (b) p(xᵢ) = r(xᵢ); ∫p = ∫qPₙ + ∫r = 0 + Σwᵢr(xᵢ) = Σwᵢp(xᵢ). (c) Π(x − xᵢ)² ≥ 0, degree 2n, positive integral, rule gives 0. (d) ℓᵢ² has degree 2n − 2, ℓᵢ(xⱼ) = δᵢⱼ, so Σwⱼℓᵢ(xⱼ)² = wᵢ = ∫ℓᵢ² > 0. (e) `x, w = gauss_legendre_nodes(5)`; `w @ x**9` (exact 0), `w @ x**10` (exact 2/11); `w.sum()`.
