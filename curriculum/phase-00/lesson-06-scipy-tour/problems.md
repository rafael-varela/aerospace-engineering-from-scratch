# P00.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and a tolerance or uncertainty on every numerical answer.

## 1. (Warm-up: three calls, three checks)

With `a = 60 N/s`, `τ = 0.3 s`, `t_end = 2 s`: compute the total impulse with `quad`, the peak with `minimize_scalar`, and the time thrust falls to 0.5 N with `brentq`. Beside each, write the analytic or hand value and the relative difference.

Answer form: three numbers with units, three reference values, three relative differences.

## 2. (Derivation: the impulse integral and its tail)

(a) Derive `∫₀ᵀ a t e^{−t/τ} dt = a τ² [1 − e^{−T/τ}(1 + T/τ)]` by parts, checking units at each line.
(b) Show that the fraction of total impulse delivered after time `T` is `e^{−T/τ}(1 + T/τ)` and find the `T/τ` at which 99 % has been delivered.
(c) Derive the composite trapezoidal rule's error bound `|E| ≤ (T h² / 12) max|F''|` from the single-interval error, and evaluate it for `n = 81` samples of this model. Compare with the error your code measures.

Answer form: derivation, one dimensionless number, one bound in N·s next to one measured error.

## 3. (Unfamiliar scenario: a different thrust shape)

A hybrid motor's thrust is modelled as a ramp to `F₀` over `t_r`, a plateau to `t_b`, then a linear tail-off to zero over `t_f`: `F₀ = 40 N`, `t_r = 0.1 s`, `t_b = 1.5 s`, `t_f = 0.2 s`. Write it as a function, compute the impulse analytically (areas of trapezia) and with `quad`, explain any disagreement, then sample it at 40 Hz and integrate with the trapezoid. Report the motor class. Finally fit `thrust_model` (the `a t e^{−t/τ}` shape) to the samples with `curve_fit` and report how badly the fitted impulse is biased.

Answer form: a function, three impulses with units, a class letter, a bias in percent with one sentence on its cause.

## 4. (Numerical, uses your code: convergence on a seed of your own)

Choose `a` and `τ` from the last two digits of your postcode or phone number (`a` = 30 + digits, in N/s; `τ` = 0.1 + digits/100, in s). Run `trapezoid_convergence` for `n = 11, 21, …, 641` and report the observed orders. Then compute `peak_thrust` and `burnout_time` to 1 N and compare each with the analytic value. State how many digits of each number you would quote in a report and why.

Answer form: the two parameters, a table of errors and orders, two comparisons, a digits decision.

## 5. (Numerical, uses your code: interpolation error budget)

The club's launch site is at 1 650 m. Compute the pressure there with all three interpolation schemes from the 2 km table and with `isa_pressure_exact`. Convert each to density using the ISA temperature (`T = 288.15 − 0.0065 h` K, `R = 287.05 J/(kg·K)`) and to a drag force `D = ½ ρ V² S C_D` on a rocket with `S = 0.0020 m²`, `C_D = 0.45` at `V = 150 m/s`. Which scheme's error in `D` is acceptable if the drag model itself is only good to 10 %? Would your answer change at 10 500 m with a 4 km table?

Answer form: a table of pressure, density and drag for four methods, two judgments.

## 6. (Judgment: which submodule, and what it will not tell you)

For each question, name the SciPy submodule and function you would reach for, the Phase 03 lesson where you will build that method, and the one thing the call will not tell you that you must check yourself:
(a) the altitude at which a rocket's velocity is zero, given a simulated trajectory table;
(b) the probability that a motor's impulse exceeds 10 N·s given a mean and standard deviation from twenty firings;
(c) the time history of altitude from a thrust curve and a drag model;
(d) the `C_D0` and `k` of a drag polar from nine wind-tunnel points;
(e) the modal frequencies of a wing idealised as a 50-node beam.

Answer form: five rows of (submodule.function, Phase 03 lesson, the unchecked assumption).

## 7. (Hard: the balance with a bad cell)

The aft cell of the balance in the lesson drifts: it reads `r_aft + 0.5 N`. Using `resolve_balance`, compute the error this produces in `L`, `D` and `M` for the true loads `[120 N, 8 N, 3 N·m]`. Then show analytically, from `B⁻¹`, that the moment error is `−0.5 ℓ` N·m regardless of the load, and that the lift error is independent of `ℓ`. Finally propose a fourth measurement that would let you *detect* a drifting cell, write the 4×3 system, and solve it in the least-squares sense with `scipy.linalg.lstsq`; show that the residual is zero for consistent readings and non-zero when the aft cell drifts.

Answer form: three numerical errors, a derivation, a 4×3 matrix, two residual norms.
