# P00.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes.

## 1. Three calls, three checks

**Assumptions.** `a = 60 N/s`, `τ = 0.3 s`, `t_end = 2 s`, the exact model (no noise).

**Result.**

| Quantity | SciPy call | SciPy value | Reference | Rel. diff. |
|---|---|---|---|---|
| Total impulse | `quad` | 5.347 312 96 N·s | `aτ²[1 − e^{−6.667}(1 + 6.667)]` = 5.347 312 96 N·s | < 10⁻¹⁰ |
| Peak | `minimize_scalar` | 6.621 83 N at 0.300 000 s | `aτ/e` = 6.621 83 N at `τ` = 0.3 s | 10⁻¹⁰ (time), < 10⁻¹² (value) |
| Time to 0.5 N | `brentq` on `[τ, t_end]` | 1.571 94 s | substitute back: `60 × 1.57194 × e^{−5.2398}` = 0.500 00 N | < 10⁻⁸ |

**Check.** The root is after the peak and before `t_end`, as the bracket requires; the impulse is below the `t → ∞` limit `aτ² = 5.40 N·s`, as it must be.

**Common mistakes.** Passing `thrust_model` (not its negative) to `minimize_scalar` and getting `t = 0`; bracketing `brentq` on `[0, t_end]`, where the function has the same sign at both ends (thrust is below 0.5 N at `t = 0` *and* at `t = 2 s`) so the call raises or finds the wrong crossing.

## 2. The impulse integral and its tail

**(a)** With `u = t`, `dv = e^{−t/τ} dt`, so `du = dt`, `v = −τ e^{−t/τ}`:

```
∫₀ᵀ t e^{−t/τ} dt = [−τ t e^{−t/τ}]₀ᵀ + τ ∫₀ᵀ e^{−t/τ} dt
                 = −τ T e^{−T/τ} + τ · τ (1 − e^{−T/τ})
                 = τ² [1 − e^{−T/τ} − (T/τ) e^{−T/τ}] = τ² [1 − e^{−T/τ}(1 + T/τ)].
```

Multiply by `a`. Units line by line: `t e^{−t/τ}` is seconds, its integral is s², times `a` in N/s gives N·s.

**(b)** Total impulse (T → ∞) is `aτ²`; the fraction delivered after `T` is `1 − [1 − e^{−x}(1 + x)] = e^{−x}(1 + x)` with `x = T/τ`. Setting it to 0.01 and solving numerically (`brentq` on `[1, 20]`): **x = 6.64**. For `τ = 0.3 s` that is 1.99 s, which is why the lesson samples to 2 s: it captures 99 % of the impulse.

**(c)** A single trapezium over `[t_i, t_i + h]` has error `−(h³/12) F''(ξ_i)`. Summing `T/h` of them: `|E| ≤ (T/h)(h³/12) max|F''| = (T h²/12) max|F''|`. Here `F'' = (a/τ) e^{−t/τ} (t/τ − 2)`, whose magnitude is largest at `t = 0`: `|F''(0)| = 2a/τ = 400 N/s²`. With `n = 81`, `h = 0.025 s`: bound `= 2 × 0.025² × 400 / 12 = 0.042 N·s`. The measured error is 0.0031 N·s, thirteen times smaller, because `|F''|` is near its maximum only in the first few hundredths of a second and is tiny over most of the range. The bound is a bound, not an estimate; both are second order in `h`.

**Check.** Dimensions of the bound: `s · s² · N/s² = N·s`. Halving `h` quarters both bound and measured error (observed rates 1.98–2.00).

**Common mistakes.** Dropping the factor `(1 + T/τ)` and reporting 5.39 N·s; using the maximum of `|F|` instead of `|F''|` in the bound; treating the bound as the expected error and concluding that 81 samples give only 1 % accuracy when they give 0.06 %.

## 3. A different thrust shape

**Assumptions.** Ramp 0–0.1 s to 40 N, plateau to 1.5 s, linear tail-off to zero at 1.7 s.

**Analytic impulse.** Triangle + rectangle + triangle: `½ × 40 × 0.1 + 40 × 1.4 + ½ × 40 × 0.2 = 2 + 56 + 4 = 62.0 N·s`. **Class F** (40–80 N·s).

**`quad`.** With `points=[0.1, 1.5, 1.7]`: 62.000 000 N·s, error estimate 7 × 10⁻¹³. Without `points` SciPy also returns 62.000 000 here: the adaptive algorithm happened to place subinterval boundaries well. Do not rely on that; with corners at less convenient times `quad` typically warns ("roundoff error is detected" or "maximum number of subdivisions") and loses digits. Passing the corners is the correct practice because the integrand is not smooth there and the Gauss–Kronrod rule inside `quad` assumes smoothness.

**Trapezoid at 40 Hz** (`h = 0.025 s`, 69 samples): 62.000 N·s exactly, because 0.1, 1.5 and 1.7 s are all multiples of 0.025 s, so every corner is a sample and every piece is linear between samples, which the trapezoid integrates exactly. Move the ramp end to 0.11 s (analytic 61.8 N·s) and the trapezoid gives 61.773 N·s, an error of −0.027 N·s = −0.04 %: the single trapezium straddling the corner is wrong by order `h × F₀ × (fraction)`. Still irrelevant beside load-cell noise.

**Fit.** `curve_fit(thrust_model, ...)` on the sampled hybrid curve converges to `a = 205 N/s`, `τ = 0.618 s` (peak 46.7 N at 0.62 s, where the real curve is a 40 N plateau). The fitted impulse to 1.7 s is 59.7 N·s, **biased low by 3.7 %**. Cause: the `a t e^{−t/τ}` family has no plateau; least squares compromises between overshooting early and undershooting late, and the area lost under the plateau exceeds the area gained at the fake peak. The class (F) survives, but the peak thrust is overstated by 17 %, which would mis-size the launch-rod loads.

**Check.** Dimensional: all impulses in N·s. Limiting case: shrink the plateau to zero (`t_b = t_r`) and the hybrid becomes a triangle whose impulse `½ F₀ (t_r + t_f)` the trapezoid still gets exactly when corners are sampled.

**Common mistakes.** Using `np.where` chains that evaluate `F₀ t / t_r` outside its range and produce negative thrust before 0; reporting the fitted peak as if it were measured; concluding that `quad` "handles corners" because one example worked.

## 4. Convergence on a seed of your own

**Example** with digits 47: `a = 77 N/s`, `τ = 0.57 s`.

| n | error [N·s] | observed order |
|---|---|---|
| 11 | 2.74 × 10⁻¹ | |
| 21 | 6.89 × 10⁻² | 1.994 |
| 41 | 1.72 × 10⁻² | 1.998 |
| 81 | 4.31 × 10⁻³ | 2.000 |
| 161 | 1.08 × 10⁻³ | 2.000 |
| 321 | 2.69 × 10⁻⁴ | 2.000 |
| 641 | 6.74 × 10⁻⁵ | 2.000 |

`peak_thrust` → 16.146 N at 0.5700 s; analytic `aτ/e = 16.146 N` at 0.57 s: agreement to 10⁻⁹ in time, 10⁻¹² in value.

`burnout_time(77, 0.57, 1.0, 2.0)` **raises `ValueError`**: at 2 s the thrust is still 4.61 N, so the threshold has not been reached and there is no root in the bracket. That is the function working correctly. With `t_end = 4 s`: 3.126 s. The impulse to 2 s is 21.6 N·s (class E) while the full impulse `aτ²` is 25.0 N·s: this slower motor delivers only 87 % of its impulse in the 2 s window, so the sampling window must be extended before the class is quoted.

**Digits.** Quote the impulse to three significant figures: on real data the noise contributes about 1 % (lesson part 5), far above the 10⁻⁵ trapezoid error, so a fourth digit is noise. The analytic peak time is exact, but on sampled data it is known only to about half a sample spacing (±0.012 s at 40 Hz), so quote 0.57 s, not 0.5700 s.

**Common mistakes.** Catching the `ValueError` and returning `t_end`; quoting digits from the printout rather than from the error analysis; forgetting that the class depends on the window when the tail is long.

## 5. Interpolation error budget

**Assumptions.** `T = 288.15 − 0.0065 h` K, `R = 287.05 J/(kg·K)`, `V = 150 m/s`, `S = 0.0020 m²`, `C_D = 0.45`. Relative error in `D` equals relative error in `p` because `ρ ∝ p` at fixed `T` and `D ∝ ρ`.

**At 1 650 m, 2 km table** (`T = 277.4 K`):

| Method | p [Pa] | ρ [kg/m³] | D [N] | error in p |
|---|---|---|---|---|
| exact | 83 011 | 1.0424 | 10.554 | 0 |
| linear | 83 315 | 1.0462 | 10.593 | +0.37 % |
| log-linear | 82 943 | 1.0416 | 10.546 | −0.08 % |
| cubic | 83 011 | 1.0424 | 10.554 | < 0.001 % |

All three are acceptable against a 10 % drag model; even linear is 27 times below the model's own uncertainty.

**At 10 500 m, 4 km table** (`T = 219.9 K`):

| Method | p [Pa] | ρ [kg/m³] | D [N] | error in p |
|---|---|---|---|---|
| exact | 24 474 | 0.3877 | 3.926 | 0 |
| linear | 25 431 | 0.4029 | 4.079 | +3.9 % |
| log-linear | 24 305 | 0.3850 | 3.899 | −0.69 % |
| cubic | 24 422 | 0.3869 | 3.917 | −0.21 % |

Linear interpolation at 3.9 % is still below 10 %, but it is now 40 % of the drag model's error budget and no longer negligible; log-linear and cubic remain negligible. The judgment changes from "any scheme" to "not linear" when the table is coarse and the function strongly curved.

**Check.** 1 650 m is near Denver's altitude; 83 kPa and 1.04 kg/m³ match the usual figures. Dimensional: `Pa / (J/(kg·K) · K) = Pa·kg/J = kg/m³`.

**Common mistakes.** Using sea-level temperature for the density conversion (a 4 % error at 1 650 m, larger than any interpolation error); comparing absolute pressure errors in pascals across altitudes where the pressures differ by a factor of three.

## 6. Which submodule, and what it will not tell you

| | Submodule.function | Phase 03 lesson | What you must check yourself |
|---|---|---|---|
| (a) | `interpolate.CubicSpline` of `v(t)` then `.roots()`, or `optimize.brentq` on the interpolant | 03 Interpolation & Splines, 02 Root Finding | That there is exactly one zero crossing in the bracket (apogee), and that the table is dense enough near it for the interpolant to be trusted. |
| (b) | `stats.norm.sf(10, loc=mean, scale=std)` | 17 Uncertainty Quantification | That twenty firings justify a normal distribution and that `std` is the population spread, not the standard error of the mean. |
| (c) | `integrate.solve_ivp` with `F(t)` and `D(v)` in the right-hand side | 05 ODE Integrators I, 06 Adaptive & Stiff | Step-size convergence (re-run with tighter `rtol` and compare), and that burnout and apogee events are detected rather than stepped over. |
| (d) | `linalg.lstsq` on `[1, C_L²]` columns (linear in the parameters) or `optimize.curve_fit` | 08 Direct Linear Solvers, 14 Unconstrained Optimization | Residual pattern: a parabola fitted to data with a stall region will show systematic residuals, and the fit will be wrong without complaint. |
| (e) | `linalg.eigh` on the stiffness and mass matrices from a beam model | 13 Finite Element Basics, 08 Direct Linear Solvers | Mesh convergence of the frequencies (50 nodes may be too few for the fourth mode) and the boundary conditions encoded in the matrices. |

**Common mistakes.** Reaching for `curve_fit` when the model is linear in its parameters (lstsq is exact and needs no initial guess); treating `solve_ivp`'s default tolerances as "accurate".

## 7. The balance with a bad cell

**Assumptions.** `ℓ = 0.15 m`; aft cell offset `δr = [0, +0.5, 0]` N; `B` as in the lesson.

**Numerical.** `resolve_balance(B x_true + δr) − x_true = [+0.500 N, −0.010 N, −0.075 N·m]` for `x_true = [120, 8, 3]`.

**Analytic.** Invert `B` by hand. From the first two rows, `L = r_fwd + r_aft` and `M = ℓ (r_fwd − r_aft)`; from the third, `D = r_drag − 0.02 L`. So

```
B⁻¹ = [  1       1      0 ]
      [ −0.02  −0.02    1 ]
      [  ℓ     −ℓ       0 ]
```

and `δx = B⁻¹ δr` is the second column times 0.5: `δL = +0.5` (independent of `ℓ` and of the load), `δD = −0.02 × 0.5 = −0.01`, `δM = −0.5 ℓ = −0.075 N·m`. The moment error is proportional to the arm because the moment is the *difference* of two readings times the arm; the lift error is the offset itself because lift is their *sum*.

**A fourth measurement.** Add a cell that reads total normal force directly at the reference point (a fourth row `[1, 0, 0]`), making `B₄` 4×3. For consistent readings `r₄ = B₄ x_true`, `lstsq` returns `x_true` with residual norm 7 × 10⁻¹⁴ (zero to rounding). With the aft drift, the residual norm is **0.289 N** and the residual vector is `[−0.167, −0.167, 0, +0.167]`: the system cannot satisfy `r_fwd + r_aft = L` and `r_total = L` at once, and the mismatch is spread over the three vertical equations. A residual above the cells' noise flags a fault; it does not, by itself, say *which* cell. (A second drag cell as the fourth row would not help: it is redundant with `D` only and leaves the vertical equations unchecked, residual 0.007 N.)

**Check.** Dimensional: row 3 of `B⁻¹` has units of length, converting N to N·m. Limiting case `ℓ → 0`: the moment error vanishes but so does the moment's observability (condition number → ∞).

**Common mistakes.** Computing `inv(B) @ r` with `np.linalg.inv` in code (fine for a 3×3 by hand, poor practice in code); choosing a redundant measurement that duplicates an existing row and then concluding that redundancy does not detect faults.

## Check-yourself answers

1. `trapezoid(y, x)`: values first, then the abscissae. Swapping them integrates `t` against `F` and returns a meaningless number with no error; the analytic comparison test is what catches it.
2. `quad` needs a callable it can evaluate anywhere and adapts its points to reach ~10⁻¹⁰; the trapezoid works on fixed samples and cannot ask for more. You cannot `quad` a sample set without first interpolating it (which adds its own error), and you would not trapezoid a model when `quad` is available and more accurate.
3. `t_peak = τ`, `F_peak = aτ/e`; `I(∞) = aτ²`.
4. `log₂(1.26 × 10⁻² / 3.15 × 10⁻³) = log₂(4.0) = 2.0`: second order. Order 1 would indicate a discontinuity in `F` inside the range (or unevenly spaced samples), and that the error falls only in proportion to `h`.
5. In the isothermal layer `ln p = ln p₁₁ − g(h − h₁₁)/(RT)` is exactly linear in `h`, so a chord through two nodes *is* the function. In the troposphere `T` varies, `ln p = ln p₀ + (g/(RL)) ln(T/T₀)` is curved, and the chord departs from it between nodes.
6. The condition number bounds the factor by which relative errors in the readings can be amplified in the recovered loads. Shrinking the arm by ten multiplies it by ten (6.7 → 67): the moment is the difference of two nearly equal readings, and the same noise gives ten times the relative moment error.
7. Raise an informative `ValueError`: no sign change means no guaranteed root in the bracket. `brentq` itself raises, but a message that says "thrust has not decayed to the threshold by `t_end`" tells the user what to change; silently widening the bracket or returning `t_end` would hide a real property of the data.
8. The optimiser found the trivial valley `a → 0` (zero thrust everywhere is a stationary point of the residual for large `τ`), usually because the initial guess was poor or absent. A data-derived guess (`τ₀` = time of the largest sample, `a₀ = F_max e / τ₀`) starts it on the right slope; bounds on `a > 0` and `τ` within the sampled range prevent the escape.
