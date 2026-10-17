# P03.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. The stability budget ($s = 15$: C_s = 1500 J/K)

**Assumptions.** Linearisation at 300 K; the radiation term contributes 4εσT³ = 4.90 W/K to the skin's loss.

**Derivation.** (a) Hand: λ_fast ≈ −(2500/1500)(1 + 1500/2 × 10⁵) = −1.6667 × 1.0075 = **−1.679 s⁻¹**; λ_slow ≈ −20/2 × 10⁵ = **−1.0 × 10⁻⁴ s⁻¹**. Code (`eigvals` of `thermal_jac`): **−1.6824 s⁻¹** and **−1.2353 × 10⁻⁴ s⁻¹**. The fast value is 0.2 % off the hand estimate (the radiation conductance); the slow one is 24 % off because the linearised radiation leaks heat in parallel with k_b (see Problem 2(c)). Stiffness ratio **1.362 × 10⁴**; slow time constant 8 100 s = 2.25 h.

(b) Explicit Euler: h < 2/1.6824 = **1.189 s**, **3 028 steps** per hour. Dormand–Prince: h < 3.3066/1.6824 = **1.965 s**, **1 832 steps** per hour (12 800 evaluations).

(c) z_fast = 60 × (−1.6824) = **−100.9**, |R| = 1/101.9 = **0.0098**: whatever skin transient exists at the start of a step is reduced to one percent by its end, so a 60 s backward-Euler step behaves like the exact solution for that mode, which has decayed to e^{−101} ≈ 0. z_slow = 60 × (−1.2353 × 10⁻⁴) = **−0.00741**, |R| = 1/1.00741 = **0.99264** against the exact e^{z} = 0.99262: the slow mode is followed with a relative error of z²/2 = 2.7 × 10⁻⁵ per step, i.e. first-order accuracy, which over 60 steps accumulates to about 0.1 K.

**Check.** Product of eigenvalues 1.6824 × 1.2353 × 10⁻⁴ = 2.078 × 10⁻⁴ s⁻²; from Problem 2(c) with radiation, det = [k_c k_b + 4εσT³(k_c + k_b)]/(C_s C_b) = (50 000 + 12 348)/(3 × 10⁸) = 2.078 × 10⁻⁴ ✓.

**Common mistakes.** Using 2/|λ| for the RK45 limit; forgetting that the step count is per hour of *simulated* time; reading |R| < 1 as "accurate" (it is only "stable").

## 2. The controller, and A- versus L-stability

(a) The estimate at step h is E = C h^q (tolerance units). At another step h_new it would be E (h_new/h)^q. Setting this to 1 gives **h_new = h E^{−1/q}**; with q = 5, h_new = h E^{−1/5}. A safety factor (0.9) is included because C is not constant (the next step starts at a different point of the solution) and because a step that lands exactly on the tolerance will be rejected half the time; aiming at 0.9 of the limit costs 0.9⁵ = 59 % of the error budget and almost no steps. The growth clamp (5) exists because the asymptotic estimate measured at h says nothing reliable about a step five times larger, and because an over-ambitious step that is rejected wastes seven evaluations.

(b) With f = λy: y_{n+1} = y_n + (h/2)(λy_n + λy_{n+1}), so y_{n+1}(1 − z/2) = y_n(1 + z/2) and

**R(z) = (1 + z/2)/(1 − z/2) = (2 + z)/(2 − z).**

|R| < 1 iff |z + 2| < |z − 2| iff z is nearer to −2 than to +2 iff Re z < 0: the trapezoidal rule is A-stable. But as z → −∞ along the real axis, R → (z/2)/(−z/2) = **−1**: a fast mode is neither damped nor amplified, only flipped in sign each step. Backward Euler has R = 1/(1 − z) → 0 (L-stable). For hλ = −300: trapezoidal |R| = 298/302 = 0.987, so a 100 K skin transient becomes −98.7 K after one step, +97.4 K after two, and survives for hundreds of steps as a sawtooth; backward Euler gives |R| = 1/301 and the transient is gone in one step. Problem 5 shows this happening.

(c) A = [[−k_c/C_s, k_c/C_s], [k_c/C_b, −(k_c + k_b)/C_b]], tr A = −k_c/C_s − (k_c + k_b)/C_b, det A = k_c(k_c + k_b)/(C_s C_b) − k_c²/(C_s C_b) = **k_c k_b/(C_s C_b)**. For the lesson's numbers tr A = −5 − 0.0126 = −5.0126 s⁻¹ and det A = 5 × 10⁴/10⁸ = 5 × 10⁻⁴ s⁻². Since 4 det/tr² = 8 × 10⁻⁵ ≪ 1, λ_fast ≈ tr A = **−5.0126 s⁻¹** and λ_slow ≈ det/tr = **−9.975 × 10⁻⁵ s⁻¹**. The computed values −5.0223 and −1.2414 × 10⁻⁴ include the radiation conductance g_r = 4εσT³ = 4 × 0.8 × 5.670 × 10⁻⁸ × 300³ = 4.90 W/K, which adds g_r/C_s = 0.0098 s⁻¹ to |λ_fast| (−5.0126 − 0.0098 = −5.0224 ✓) and changes the determinant to [k_c k_b + g_r(k_c + k_b)]/(C_s C_b) = (50 000 + 12 348)/10⁸ = 6.235 × 10⁻⁴, so λ_slow ≈ −6.235 × 10⁻⁴/5.0224 = −1.2414 × 10⁻⁴ ✓.

**Common mistakes.** Writing R_trap = 1 + z + z²/2 (that is its Taylor expansion, not R); claiming A-stability implies damping of fast modes; dropping the k_c² term from the determinant.

## 3. Actuator in the loop ($s = 15$: τ = 15 ms)

**Derivation.** (a) λ_act = −1/τ = **−66.7 s⁻¹**; λ_phugoid = ±2πi/60 = **±0.1047i s⁻¹**; λ_sp = ±2πi/3 = **±2.094i s⁻¹**. Stiffness ratio (by |λ|, since the aircraft modes are oscillatory) 66.7/0.1047 = **637**.

(b) h_stab = 3.307 × 0.015 = **0.0496 s**; 600/0.0496 = **12 100 steps**, about **84 700 evaluations** (72 600 with FSAL). The accuracy step for the short-period mode at rtol = 10⁻⁶ is about 3 s/40 = 0.075 s, i.e. 8 000 steps: for τ = 15 ms the problem is only mildly stiff (the stable step is two thirds of the accuracy step), and RK45 would cost about 1.5× what accuracy alone would need. For τ = 5 ms the stable step falls to 0.0165 s and the ratio to 4.5: that is where it starts to hurt.

(c) With the lag removed, the eigenvalues are ±0.1047i and ±2.094i: ratio 20, not stiff. Lost: the actuator's phase lag at the short-period frequency (φ = atan(ωτ) = atan(2.094 × 0.015) = 1.8° at τ = 15 ms, small; 10° at τ = 85 ms) and its rate limit, which is a nonlinear effect no eigenvalue shows and which decides whether the control law saturates.

(d) `Radau` (or `BDF`), supplying `jac`, since Newton with a finite-difference Jacobian of an n-state flight model costs n extra evaluations per Jacobian; also a `first_step` of the order of τ so that the actuator's initial transient is resolved rather than stepped over.

**Check.** Dimensionally, τ in s gives λ in s⁻¹; a period of 3 s and 40 steps per period give 0.075 s ✓.

**Common mistakes.** Using the real parts only (zero for the undamped modes, giving an infinite ratio that says nothing); forgetting that stiffness depends on the integration *length*: over 0.1 s of flight nothing here is stiff.

## 4. Your own entry angle ($s = 15$: γ_E = −5.5°)

**Assumptions.** The lesson's capsule (β = 222.2 kg/m²), exponential atmosphere, entry at 120 km and 7 800 m/s.

**Derivation.** (a) |sin 5.5°| = 0.09585; a_max = 7800² × 0.09585/(2e × 7200) = 5.832 × 10⁶/39 143 = 149.0 m/s² = **15.19 g**; h* = 7200 ln[1.225 × 7200/(222.2 × 0.09585)] = 7200 ln 414.1 = **43.39 km**.

(b) With gravity (`atol = 1e-3`):

| rtol | accepted | rejected | nfev | peak read (g) | at (km) | t (s) |
|---|---|---|---|---|---|---|
| 10⁻⁶ | 28 | 9 | 261 | 16.81 | 42.15 | 107.7 |
| 10⁻⁸ | 42 | 9 | 359 | 16.72 | 41.38 | 109.2 |

(c) The finely resolved run (rtol 10⁻¹⁰, h_max = 0.25 s) gives **16.81 g at 42.16 km, t = 107.6 s**. The 10⁻⁶ run happened to place a step within 0.1 s of the peak; the 10⁻⁸ run's nearest stored point is 1.6 s late and reads 0.09 g low. The *integration* at 10⁻⁸ is a hundred times more accurate; the *reading* is worse because the peak is sampled, not located. Gravity lifts the peak from 15.19 g (drag only) to 16.81 g and lowers it from 43.4 to 42.2 km.

(d) log(42/28)/log 100 = **0.088**, against the theoretical 0.2. Of the 28 steps at 10⁻⁶, five are start-up steps (0.02 → 0.1 → 0.5 → 2.4 → 12 s, growing by the clamp of 5 each time) and nine attempts were rejected; those counts are nearly the same at both tolerances, so only the error-limited steps (23 against 37) scale with rtol, and their ratio 1.6 gives an exponent of 0.10. The remainder is the clamp and the start-up: the −1/5 law assumes every step is chosen by the error estimate alone.

**Check.** a_max/|sin γ_E| is the same for −5.5° and −6° (15.19/0.09585 = 158.5 = 16.57/0.1045 ✓).

**Common mistakes.** Reporting the stored-point peak as "the" peak without a resolution statement; forgetting to convert m/s² to g with 9.80665.

## 5. A second-order implicit method

(a) With f = λy and z = hλ: G(z_new) = z_new − y − (h/2)(λy + λz_new) = 0 gives z_new(1 − z/2) = y(1 + z/2), so one step returns y(1 + hλ/2)/(1 − hλ/2) ✓, and Newton reaches it in one iteration because G is linear.

(b) From T_s(0) = 400 K, T_b(0) = 300 K; Radau reference after one hour: T_s = 329.653 K, T_b = 328.922 K.

| h (s) | trapezoidal T_s at 0, h, 2h, 3h (K) | final error T_s, T_b (K) | backward Euler T_s at 0, h, 2h, 3h (K) | final error T_s, T_b (K) |
|---|---|---|---|---|
| 60 | 400, 203.7, 398.4, 207.3 | 42.9, 1.89 | 400, 302.0, 302.2, 302.8 | 0.086, 0.086 |
| 30 | 400, 204.7, 395.4, 210.3 | 3.62, 0.49 | 400, 302.0, 301.6, 301.9 | 0.043, 0.043 |
| 15 | 400, 207.2, 390.3, 216.9 | 0.11, 0.11 | 400, 302.5, 301.4, 301.5 | 0.022, 0.022 |

(c) The 100 K excess of the skin is flipped, not damped: with z = 60 × (−5.022) = −301.3, R_trap = (2 + z)/(2 − z) = **−0.987**, so the 100 K excess becomes −98.7 K after one step (T_s ≈ 201 K, or 203.7 K once the structure's small response is included) and +97.4 K after two (398.4 K). The sawtooth decays by only 1.3 % per step at h = 60 s and survives the whole hour (42.9 K error at the end); at 15 s, |R| = 0.949 and it has decayed by 0.949²⁴⁰ ≈ 3 × 10⁻⁶ by the end, so the final error is small. The final errors (42.9, 3.62, 0.11 K) fall far faster than h²: they are the surviving transient, not the truncation error, and no order can be read from them. Backward Euler, with |R| = 1/302, removes the excess in one step and shows its first-order error (0.086, 0.043, 0.022 K: ∝ h).

(d) At h = 1, 0.5, 0.25 s the trapezoidal errors are 5.42 × 10⁻⁴, 1.52 × 10⁻⁴, 4.10 × 10⁻⁵ K: observed order **1.86**, approaching 2. The asymptotic order appears once |hλ_fast| is of order 1 or less (h ≲ 0.2 s), that is, once the fast mode is actually resolved; a non-L-stable method must resolve the transient to earn its order, which is the stability limit of an explicit method in disguise.

**Common mistakes.** Starting Newton from the explicit predictor (from 400 K at h = 60 s it predicts T_s = −29 500 K and the T⁴ term makes Newton diverge); concluding from (b) that "the trapezoidal rule is unstable" (it is stable: nothing grows).

## 6. Choosing the integrator for the coupled simulation

(a) The thermal block: the skin eigenvalue −5.02 s⁻¹ is unchanged by the coupling (the trajectory does not depend on the temperatures, so the Jacobian is block-triangular and the eigenvalues are those of the blocks). The trajectory's own rates are of order v̇/v ≈ 0.04 s⁻¹ at the peak. The stable RK45 step is **0.66 s** for the whole coupled system.

(b) (i) 200/0.66 ≈ 300 steps × 7 ≈ **2 100 evaluations** per case (plus start-up); (ii) 200 steps × ≈ 3 (Newton with an analytic Jacobian) ≈ **600 evaluations** plus 200 Jacobians and 6×6 solves; (iii) Radau: the lesson's 17 steps per smooth hour will become perhaps 100 steps through the heating pulse, ≈ **800 evaluations**. At 20 μs per evaluation and 10 000 cases: 420 s, 120 s and 160 s. All three fit an overnight budget by a wide margin; the budget does not discriminate, so the choice is on robustness and accuracy: Radau (fifth order, L-stable, error-controlled) over backward Euler (first order), with RK45 as the cross-check.

(c) `atol` is an array: 10⁻³ m and 10⁻⁶ m/s on the trajectory states with rtol = 10⁻⁸, and 0.01 K on the two temperatures. The error norm (1) then gives each component its own meaning of "small enough". Alternatively, because the coupling is one-way, integrate the trajectory first at 10⁻⁸ and feed q(t) through an interpolant to a separate thermal integration, which also allows the two solvers of (d) without any splitting error, since the trajectory is known in full before the thermal solve starts.

(d) Gain: each solver keeps its natural step (RK45 about 10 s on the trajectory away from the peak, backward Euler 1 s on the thermal model), and the trajectory code is untouched. Risk: a one-way exchange every 1 s is a first-order splitting error in q, and near the peak q changes by 20 % per second; the skin temperature could be off by several kelvins. Test: halve the exchange interval to 0.5 s and 0.25 s; if the peak skin temperature moves by less than the 20 % input uncertainty implies (tens of kelvins), the splitting is harmless; if it moves comparably, interpolate q within the interval or use the one-way ordering of (c).

## Check-yourself answers

1. The fifth-order solution is more accurate, so it is used to advance (local extrapolation); the estimate y₅ − y₄ is the error of y₄ and is a conservative bound on the error actually committed.
2. The estimate is O(h⁵) (local error of the fourth-order solution), so err ∝ h⁵ and the step giving a target error scales as (target/err)^{1/5}; a 2(3) pair has an O(h³) estimate and uses −1/3.
3. Explicit Euler R = 1 + z, interval (−2, 0); RK4 R = 1 + z + z²/2 + z³/6 + z⁴/24, interval (−2.785, 0); implicit Euler R = 1/(1 − z), stable on the whole negative axis (and the whole left half-plane).
4. 4.05 × 10⁴; h < 2/5.022 = 0.398 s; at least 9 050 steps per hour (the lesson's run shows 0.3 s is stable, 0.5 s overflows).
5. a_max = v_E²|sin γ_E|/(2eH), independent of β; h* = H ln[ρ₀H/(β|sin γ_E|)], which falls as β rises (a denser vehicle goes deeper before peaking).
6. The accepted step stops responding to `rtol`: tightening the tolerance by 100 changes the step count by far less than 2.5×, and the error delivered is orders of magnitude below the tolerance (10⁻⁷ K at rtol 10⁻⁶ in the lesson). A rejection pattern of "accept, accept, reject" at a fixed step size is the same symptom seen close up.
7. The explicit predictor y_n + h f(t_n, y_n) moves by h times the fast rate, which for h ≫ 1/|λ_fast| is a huge overshoot (−29 500 K in Problem 5); starting from y_n is within the fast mode's own basin, and on a linear problem Newton converges in one step from anywhere.
8. The relative scale rtol·|y| is zero for a component that is zero, so the controller would demand infinite accuracy there; `atol` sets the error below which a component is "small enough". With atol = 10⁻⁹ m on the range coordinate (which starts at 0) the first step was 0.2 μs; with 10⁻³ m it was 0.02 s.
