# P00.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `integrate.quad(lambda t: thrust_model(t, a, tau), 0, t_end)[0]` against `analytic_impulse`. 2. `optimize.minimize_scalar(lambda t: -thrust_model(...), bounds=(0, t_end), method="bounded")` against `(τ, aτ/e)`. 3. `optimize.brentq(lambda t: thrust_model(...) - 0.5, tau, t_end)`; check by substituting back into the model. 4. Report `|num − ref| / ref` for each.

## 2.
1. (a) Parts with `u = t`, `dv = e^{−t/τ} dt`, `v = −τ e^{−t/τ}`; evaluate the boundary term and the remaining integral; collect to the bracket form; units `(N/s)(s²) = N·s`. 2. (b) Tail fraction = `1 − bracket`; solve `e^{−x}(1 + x) = 0.01` numerically (`brentq`). 3. (c) Sum `T/h` single-interval errors; bound `|F''|` by evaluating `F'' = a e^{−t/τ}(t/τ − 2)/τ` and taking its maximum magnitude on `[0, T]`, which is at `t = 0`; insert `h = T/(n−1)`. Compare with `trapezoid_convergence` at `n = 81`.

## 3.
1. Define the piecewise function with `np.piecewise` or nested `np.where`. 2. Analytic impulse: triangle `½ F₀ t_r` + rectangle `F₀ (t_b − t_r)` + triangle `½ F₀ t_f`. 3. `quad` with `points=[t_r, t_b, t_b + t_f]`; without `points` note any warning and the error. 4. Sample at 40 Hz from 0 to `t_b + t_f`, trapezoid; the samples straddle the corners, so expect an error of order `h × F₀`. 5. Class from `motor_class`. 6. `curve_fit(thrust_model, t, f, p0=...)`; compare `analytic_impulse(a_fit, tau_fit, t_end)` with the true impulse; the bias comes from fitting a peaked shape to a plateau.

## 4.
1. Derive `a`, `τ` from your digits and state them. 2. `trapezoid_convergence(a, tau, 2.0, [11, 21, 41, 81, 161, 321, 641])`; orders from successive ratios. 3. `peak_thrust` vs `(τ, aτ/e)`; `burnout_time(..., 1.0, 2.0)` vs `brentq` on the exact model (they should agree to the solver tolerance). 4. Digits: the impulse from 81 noisy samples is good to about 1 %, so three significant figures; the analytic peak is exact, but a real peak would be limited by sample spacing and noise.

## 5.
1. `interpolate_pressure(np.array([1650.0]), kind)` for the three kinds, plus `isa_pressure_exact`. 2. `ρ = p / (R T)` with `T` from the lapse rate. 3. `D = ½ ρ V² S C_D`. 4. The relative error in `D` equals that in `p`; compare each with 10 %. 5. Repeat at 10 500 m with a 4 km table (build `hn = np.arange(0, 20001, 4000)` and interpolate as in the lesson) and compare again.

## 6.
1. For each: identify whether the task is a root, an initial-value ODE, a probability from a distribution, a linear or nonlinear least-squares fit, or an eigenvalue problem. 2. Map to `optimize.brentq` (or `interpolate` + root), `integrate.solve_ivp`, `stats.norm`, `linalg.lstsq` (or `optimize.curve_fit`), `linalg.eigh`. 3. Name the Phase 03 lesson (02 Root Finding, 05/06 ODE Integrators, 17 Uncertainty Quantification, 08 Direct Linear Solvers or 14 Unconstrained Optimization, 08/09 solvers for eigenproblems) and the assumption the call takes on faith.

## 7.
1. `resolve_balance(B @ x_true + [0, 0.5, 0]) − x_true`. 2. Compute `B⁻¹` symbolically from the structure of `B` (two rows share `L/2`, differ in sign of the moment term); `δx = B⁻¹ δr` picks out the second column of `B⁻¹`. 3. Add a fourth cell, e.g. a second drag cell or a vertical cell at the reference point, giving a 4×3 `B₄`. 4. `linalg.lstsq(B₄, r)`; residual `‖B₄ x − r‖` is zero for `r = B₄ x_true` and non-zero with the drift; the size of the residual flags the fault.
