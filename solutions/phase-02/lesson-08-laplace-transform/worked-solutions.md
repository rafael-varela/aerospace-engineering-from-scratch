# P02.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. All simulated records are synthetic.

## 1. A transform and an inverse by hand ($s = 15$: a = 3)

**Derivation and result.** (a) By linearity and the table, F(s) = 3/s + 2/(s + 3) + 4/(s² + 16). At s = 2: 1.5 + 0.4 + 4/20 = **2.1**; `laplace_numeric` returns 2.1000000. (b) s² + 4s + 13 = (s + 2)² + 9, so Y = (s + 2)/((s + 2)² + 3²), the damped-cosine entry with a = 2, ω = 3: **y(t) = e^{−2t} cos 3t**. `residues([1, 2], [1, 4, 13])`: poles **−2 ± 3i**, residues **0.5 and 0.5**; so y = 0.5e^{(−2+3i)t} + 0.5e^{(−2−3i)t} = e^{−2t} cos 3t, the same. Values at t = 0, 0.3, 1.0: **1.0000, 0.3411, −0.1340**. (c) At s = 0 the integrand of L{3} is the constant 3 and ∫₀^∞ 3 dt diverges; the transform exists only for Re s > 0, and F(s) → ∞ as s → 0 through the 3/s term, which is the pole of the step.

**Check.** y(0) = 1 equals lim_{s→∞} sY(s) = 1 (initial-value theorem) ✓. Dimensions: F in (units of f) × s; 1/s has units s ✓.

**Common mistakes.** Writing 1/(s + 2)² + 9 and reading an exponential times t; forgetting the conjugate residue and reporting a complex y.

## 2. The derivative rule, the step response, and the overshoot

**Derivation.** (a) L{f′} = ∫₀^∞ f′ e^{−st} dt = [f e^{−st}]₀^∞ + s ∫₀^∞ f e^{−st} dt. If |f(t)| ≤ M e^{at} and Re s > a, f e^{−st} → 0 at ∞, leaving −f(0) + sF(s). Apply to f′: L{f″} = s L{f′} − f′(0) = s²F − s f(0) − f′(0).

(b) (s²Y − s y₀ − ẏ₀) + 2ζω_n(sY − y₀) + ω_n² Y = ω_n² U, so Y = [ω_n² U + (s + 2ζω_n) y₀ + ẏ₀]/(s² + 2ζω_n s + ω_n²).

(c) Y = ω_n²/[s(s − p)(s − p̄)] with p = −ζω_n + iω_d. Residue at 0: ω_n²/|p|² = 1. Residue at p: r = ω_n²/[p(p − p̄)] = ω_n²/[p · 2iω_d]. Then y = 1 + r e^{pt} + r̄ e^{p̄t} = 1 + 2 Re(r e^{pt}). Compute r: with p = ω_n e^{iφ} where cos φ = −ζ, sin φ = √(1 − ζ²), r = ω_n²/(2iω_n√(1 − ζ²) · ω_n e^{iφ}) = e^{−iφ}/(2i√(1 − ζ²)) = −i e^{−iφ}/(2√(1 − ζ²)). So 2 Re(r e^{pt}) = e^{−ζω_n t} Re(−i e^{i(ω_d t − φ)})/√(1 − ζ²) = e^{−ζω_n t} sin(ω_d t − φ)/√(1 − ζ²). Expand sin(ω_d t − φ) = sin ω_d t cos φ − cos ω_d t sin φ = −ζ sin ω_d t − √(1 − ζ²) cos ω_d t; dividing by √(1 − ζ²): y = 1 − e^{−ζω_n t}[cos ω_d t + (ζ/√(1 − ζ²)) sin ω_d t]. For the actuator r = −0.5 − 0.375i, as the code prints.

(d) ẏ = e^{−ζω_n t}{ζω_n[cos ω_d t + (ζ/√(1 − ζ²)) sin ω_d t] + ω_d sin ω_d t − ζω_n cos ω_d t... } collect: the cos terms cancel and the sin terms give ω_n sin ω_d t/√(1 − ζ²) [ζ² + (1 − ζ²)] = (ω_n/√(1 − ζ²)) e^{−ζω_n t} sin ω_d t. First zero after t = 0 at ω_d t = π: **t_p = π/ω_d**. There cos ω_d t = −1, sin = 0, so y(t_p) = 1 + e^{−ζω_n π/ω_d} = 1 + e^{−πζ/√(1 − ζ²)}: **M_p = e^{−πζ/√(1 − ζ²)}**. Actuator: t_p = π/20 = **157.1 ms**, M_p = e^{−0.6π/0.8} = **9.48 %**.

**Check.** ζ → 0: M_p → 1 (the undamped response reaches 2); ζ → 1: M_p → 0 and t_p → ∞ (no overshoot). y(0) = 1 − 1 = 0 and ẏ(0) = 0 ✓ (from rest).

**Common mistakes.** Dropping the −f(0) term in (a), which puts the initial conditions in the wrong place in (b); in (c), forgetting that the residue at p is complex and reporting a non-decaying answer.

## 3. Pressure lag in a pitot-static line ($s = 15$: τ = 1.25 s, climb 8.75 m/s, R = −105 Pa/s)

**Assumptions.** Linear first-order lag; constant climb rate; constant pressure gradient 12 Pa/m (true near sea level in the ISA, where dp/dh = −ρg = −12.0 Pa/m).

**Derivation.** (a) τ(sP − 0) + P = R/s², so **P_ind(s) = R/(s²(τs + 1))**: a **double pole at s = 0** and a simple pole at s = −1/τ = −0.8 s⁻¹. `residues` raises `ValueError("repeated pole …")` because N(p)/D′(p) divides by D′(0) = 0. (b) Partial fractions: R/(s²(τs + 1)) = A/s² + B/s + C/(τs + 1). Cover s²: A = R. Cover (τs + 1) at s = −1/τ: C = R/(1/τ²) = Rτ². Coefficient of s¹ in the numerator R = A(τs + 1) + Bs(τs + 1) + Cs²: Aτ + B = 0 so B = −Rτ. Inverting: **p_ind(t) = R[t − τ + τ e^{−t/τ}]**. For t ≫ τ the exponential dies and p_ind ≈ R(t − τ) = p_true − Rτ: a constant lag of **Rτ = −131.25 Pa**. SymPy: `inverse_laplace_transform` gives −105t + 131.25 − 131.25e^{−0.8t} ✓. (c) The indicated pressure is 131.25 Pa *higher* than the true pressure (R < 0), which the altimeter reads as 131.25/12 = **10.9 m lower** than the true altitude: the altimeter reads low during the climb, by the climb rate times τ (8.75 × 1.25 = 10.9 m). (d) `solve_ivp` at t = 1, 5, 10 s: **−32.72, −396.15, −918.79 Pa**; the formula gives −32.72, −396.15, −918.79 ✓. (e) E(s) = R/s² − R/(s²(τs + 1)) = R[(τs + 1) − 1]/(s²(τs + 1)) = Rτ/(s(τs + 1)); s E(s) → Rτ as s → 0, so the error tends to the constant Rτ = −131.25 Pa (the theorem applies because sE has its only pole at −1/τ, in the left half-plane).

**Check.** Dimensions: Rτ is Pa/s × s = Pa ✓. Limiting case τ → 0: no lag. The lag altitude error equals climb rate × τ, independent of the pressure gradient, as it must (the gradient converts both ways).

**Common mistakes.** Using the simple-pole formula at the double pole (the code refuses on purpose); getting the sign of the altitude error backwards (pressure falls with altitude).

## 4. Does a softer actuator meet the specification? ($s = 15$: ζ = 0.45)

**Result.** (a) Poles **−11.25 ± 22.33i s⁻¹**; overshoot **20.53 %**; peak time **140.7 ms**; 2 % settling time **333.8 ms** on the 10 µs grid, against 4/(ζω_n) = 355.6 ms (the rule bounds the envelope; this response leaves the band for the last time 6 % earlier). (b) **Fails both halves:** 20.5 % > 10 % and 334 ms > 250 ms. (c) 50 ms pulse: peak **0.671** of the command at **76.8 ms** (the pulse ends at 50 ms, before the surface has arrived). (d) |G| and lag: 3 rad/s: **1.009, 6.3°**; 10 rad/s: **1.094, 23.2°**; 25 rad/s: **1.111, 90°**. At the short-period frequency the actuator passes commands with 0.9 % gain error and 6° lag: effectively undistorted, which is why actuators are designed with ω_n several times the airframe's modes.

**Check.** Overshoot from the sampled peak: 1.2053 − 1 ✓ (test tolerance 10⁻⁴). The poles' modulus is √(11.25² + 22.33²) = 25.0 = ω_n ✓.

**Common mistakes.** Quoting 4/(ζω_n) as *the* settling time; reading the pulse peak before the pulse ends.

## 5. What the transfer function cannot tell you

(a) The 1° and 3° responses scale (9.5 % and 9.4 % overshoot) so the actuator is linear up to about 3°. At 10° the response is a 45°/s ramp with 2 % overshoot: the valve is passing its maximum flow and the surface moves at a fixed **rate limit of 45°/s**. Superposition fails (the 10° response is not ten times the 1° one), so no transfer function, which by construction obeys superposition, describes it.

(b) At 15° the rate limit dominates: the time to reach the command is about 15/45 = **0.33 s** regardless of ω_n, so peak time, settling time and bandwidth are set by the limit, not by the poles, and overshoot is small because the surface arrives slowly. Measure the maximum surface rate and the time to 90 % at full amplitude; keep the linear numbers for small signals.

(c) Overshoot < 5 % needs ζ > 0.69 (solve e^{−πζ/√(1−ζ²)} = 0.05), a wedge of half-angle arccos 0.69 = **46°** about the negative real axis; settling within 150 ms needs ζω_n > 4/0.15 = **26.7 s⁻¹**, a vertical line at Re s = −26.7. The current poles (Re s = −15) are outside. Changing ζ alone moves the poles on the circle of radius 25; even at ζ = 1 the real part is only −25 > −26.7, so **ζ alone cannot do it**; ω_n must rise to at least 26.7/0.69 = 39 rad/s.

(d) ω_n² ∝ stiffness, so ω_n falls to √0.7 × 25 = **20.9 rad/s** at constant ζ (if the damping coefficient c is unchanged, ζ = c/(2√(km)) *rises* to 0.6/√0.7 = 0.72, which helps overshoot). The poles move toward the origin; the settling time rises to 4/(0.6 × 20.9) = 319 ms at fixed ζ: the **250 ms settling item is at risk**, the overshoot is not.

**Check.** The linear step response of the lesson's actuator has a maximum rate of 0.499 ω_n per unit step (from the sampled response), so a 10° step would demand 125°/s, far above the 45°/s limit; a 1° step demands 12.5°/s and stays linear, consistent with the bench data.

**Common mistakes.** Calling a rate-limited actuator "more damped"; forgetting that ζ depends on k as well as c.

## 6. The shock-response curve of a half-sine pulse ($s = 15$: t_p = 25 ms)

**Derivation.** (a) Let f(t) = sin ω_p t for t ≥ 0 with ω_p = π/t_p. Then f(t − t_p) u_{t_p}(t) = sin(ω_p t − π) u_{t_p} = −sin(ω_p t) u_{t_p}, so f + f(t − t_p)u_{t_p} = sin ω_p t for t < t_p and 0 for t > t_p: the half-sine pulse. By the shift rule, U(s) = F(s)(1 + e^{−t_p s}) with F = ω_p/(s² + ω_p²). (b) By linearity and the shift rule, y = y_sin(t) + y_sin(t − t_p) u_{t_p}(t), where Y_sin = G · ω_p/(s² + ω_p²) has four simple poles (−15 ± 20i, ±iω_p) and `inverse_laplace` returns it as a real sum; the shifted copy is `inverse_laplace(..., t - t_p)`, zero before t_p by causality. For t_p = 25 ms (ω_p = 125.7 rad/s): **peak y = 0.1967 at t = 59.3 ms**. (c) `solve_ivp` with the pulse as input (max_step 0.5 ms): maximum difference **4.7 × 10⁻¹¹**; `convolve_causal` of the impulse response with the sampled pulse: **6.2 × 10⁻⁷** (first-order quadrature at 50 µs). (d) Peaks for ζ = 0.6:

| t_p [ms] | ω_p/ω_n | peak y | time of peak [ms] |
|---|---|---|---|
| 25 | 5.03 | 0.197 | 59 |
| 50 | 2.51 | 0.383 | 73 |
| 100 | 1.26 | 0.686 | 104 |
| 200 | 0.63 | 0.960 | 161 |
| 400 | 0.31 | 1.028 | 252 |

Short pulses: the system cannot follow and sees only the impulse ∫u dt = 2t_p/π; the peak is that area times max g = 0.0159 × 12.5 = 0.198 at t_p = 25 ms, within 1 % of 0.197, and the peak occurs at the impulse response's own peak time (59 ms). Long pulses: the input is slow compared with 1/ω_n, the surface tracks it, and the peak tends to the pulse amplitude 1 with a small overshoot that approaches the step overshoot's fraction as t_p → ∞ (1.028 at 400 ms).

**Check.** Dimensions: ω_p t_p = π ✓. The convolution and the residue solution are two different theorems (convolution rule, shift rule) giving the same curve to 10⁻⁶: a stronger check than either alone.

**Common mistakes.** Transforming the half-sine as if it were a full sine (missing the e^{−t_p s} term); evaluating `inverse_laplace` at t − t_p without the causal cut (it is built in here, but a hand version must include u_{t_p}).

## Check-yourself answers

1. F(s) = ∫₀^∞ f(t)e^{−st} dt; L{e^{−3t}} = ∫₀^∞ e^{−(s+3)t} dt = 1/(s + 3) for Re s > −3.
2. L{f′} = sF − f(0); L{f″} = s²F − s f(0) − f′(0); the initial conditions become polynomial terms added to the input side, with the same denominator as the forced response.
3. L{t} = 1/s², and the s-shift rule gives L{t e^{−3t}} = 1/(s + 3)².
4. Poles −2 ± 5i: stable (real part −2 < 0); damped frequency 5 rad/s, decay rate 2 s⁻¹ (time constant 0.5 s).
5. M_p = e^{−0.5π/√0.75} = 16.3 %; it depends on ζ only, not on ω_n.
6. t_s ≈ 4/(ζω_n) = 4/15 = 267 ms; it bounds the envelope e^{−ζω_n t}, and the actual oscillation may already be inside the band when the envelope crosses it (238 ms here).
7. e^{−2s}/s; the response is y_step(t − 2) u_2(t): the ordinary step response, delayed by 2 s and zero before.
8. Amplitude ratio 0.833, phase lag 90°: the complex gain H(r) = 1/(1 − r² + 2iζr) of P02.07 at r = 1.
