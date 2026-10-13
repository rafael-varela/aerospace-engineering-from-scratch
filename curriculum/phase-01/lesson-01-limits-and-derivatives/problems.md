# P01.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, give angles in both radians and degrees where a derivative with respect to angle appears, and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: slope of your own lift curve)

Your section has C_Lmax = 1.2 + s/50 and the thin-airfoil slope a₀ = 2π rad⁻¹, with the lesson's model C_L = C_Lmax tanh(a₀α/C_Lmax). At α = s/2 degrees, find by hand: (a) u = a₀α/C_Lmax; (b) C_L; (c) dC_L/dα in rad⁻¹ and in deg⁻¹; (d) the fraction of the thin-airfoil slope that remains. Check (c) with `central_difference(lambda x: lift_coefficient(x, cl_max=...), alpha, 1e-4)` and report the difference. Answer form: four numbers with units, four significant figures, plus the check.

## 2. (Derivation: the error of each formula)

Starting from the expansion f(x ± h) = f ± h f′ + (h²/2) f″ ± (h³/6) f‴ + (h⁴/24) f⁽⁴⁾ ± (h⁵/120) f⁽⁵⁾ + …:

(a) Derive the leading error term of the forward difference and of the central difference, and state the order of each.
(b) Write the central difference with step 2h, combine it with the step-h version so that the h² terms cancel, and show the result is the five-point stencil. What is its leading error term?
(c) Prove that Richardson extrapolation [4D(h/2) − D(h)]/3 of the central difference D is identical to the five-point stencil with step h/2.
(d) For which polynomial degree is each of the three formulas exact, and why?

Answer form: derivation with each cancellation shown; one sentence for (d).

## 3. (Unfamiliar scenario: Mach number in a climb)

An aircraft climbs through the ISA troposphere at constant true airspeed V = 150 m/s and constant climb rate ḣ = 10 m/s. The Mach number is M = V/a with a = √(γRT), γ = 1.4, R = 287.05287 J/(kg K), and T(h) from `isa_temperature_two_layer`.

(a) Use the chain rule to write dM/dt in terms of M, T, dT/dh and ḣ, without evaluating anything.
(b) Evaluate dM/dt at h = 5000 m and at h = 10 000 m. Does M rise or fall in the climb, and why physically, given that V is constant?
(c) Check (b) numerically: form M(h) as a Python function and use `central_difference` with h = 1 m, times ḣ.
(d) What happens to dM/dt as the aircraft crosses 11 km? Use `one_sided_derivatives` on M(h) at 11 000 m and explain the result in one sentence.

Answer form: one formula, four numbers with units, one sentence each for (b) and (d).

## 4. (Numerical, uses your code: the error curve for two functions)

Run `error_vs_step` for h = 10⁻¹ … 10⁻¹² on (i) sin x at x = 1 and (ii) the lesson's lift model at α = 16°.

(a) Tabulate and plot both on log-log axes (one figure, two panels). Fit `observed_order` for forward, central and five-point in the truncation regime of each, stating which h you included.
(b) Predict h_opt for the central difference in both cases with `optimal_step_central` (get f‴ from SymPy) and compare with the measured minimum.
(c) Why is the 16° curve's minimum at a smaller h than the sin curve's? Use the two h_opt values to answer.
(d) ★ Implement the complex-step derivative Im f(x + ih)/h for the lift model (`np.tanh` accepts complex arguments) and add it to the 16° plot for h = 10⁻⁴, 10⁻⁸, 10⁻¹², 10⁻²⁰. Explain in two sentences why it has no rounding floor.

Answer form: two tables, one figure, six fitted orders, two h_opt comparisons, one paragraph.

## 5. (Derivative from a table with noise)

A wind-tunnel table of C_L against α at 1° spacing from 0° to 14° is simulated by the lesson's model plus uniform noise of ±0.005 (a realistic balance resolution):

```python
rng = np.random.default_rng(s)
deg = np.arange(0, 15, 1.0)
cl_tab = lift_coefficient(np.radians(deg)) + rng.uniform(-0.005, 0.005, size=deg.size)
```

(a) Estimate dC_L/dα at 8° with central differences using the neighbours 1°, 2° and 4° away. Compare each with the exact 4.3427 rad⁻¹.
(b) For each spacing compute the predicted truncation error (h²/6)|C_L‴| with C_L‴ = −12.85 rad⁻³ and the worst-case noise error 0.005/h. Which dominates?
(c) Balance the two to find the spacing that minimises the worst-case total error, in degrees. Comment on how different it is from the 10⁻⁶ rad optimum of the smooth case.
(d) Repeat (a) with the noise removed. Which spacing is best now?

Answer form: a table of estimates and errors, one spacing in degrees, two sentences.

## 6. (Judgment: when a numerical derivative lies)

For each situation, say what goes wrong with a central difference, give one measurable sign in the output that it has gone wrong, and say what you would do instead:

(a) The function has a kink at the point (ISA temperature at 11 km).
(b) The function is a table with three significant figures.
(c) The point is just below a sharp stall (a real section, not the tanh model).
(d) The function is pressure across a shock wave in supersonic flow (a jump, not a kink).
(e) The step is 10⁻¹² rad.

Answer form: five short paragraphs, each with one number or one concrete observable.

## 7. (Hard: the inverse problem, by hand and by code)

A gust analysis asks the question backwards: how much angle of attack does it take to produce a given lift change?

(a) Using the inverse function rule, write dα/dC_L in terms of dC_L/dα, and evaluate it at the angle where C_L = 1.0 for the lesson's model (C_Lmax = 1.4). Give the result in rad and in degrees per unit C_L.
(b) Find that angle two ways: exactly, by inverting the tanh; and numerically, by using the intermediate value theorem to bracket it between 8° and 12° and bisecting to 10⁻⁶ rad. How many bisection steps are needed?
(c) Use (a) to estimate the Δα needed for ΔC_L = +0.1 from C_L = 1.0, then compute the exact Δα by inverting the model at C_L = 1.1. How large is the linearisation error, and why is it larger here than it would be at C_L = 0.5?
(d) Write one sentence on why a certification engineer would want (a) for a *measured* lift curve, and what step size they should use for it.

Answer form: one formula, four numbers with units, an integer, two sentences.
