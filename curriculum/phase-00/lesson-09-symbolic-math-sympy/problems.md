# P00.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (symbolic) and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the package)

An instrument package has mass $m = 0.5 + 0.1s$ kg, drag coefficient $C_D = 1.0$ and frontal area $A = 0.010$ m². Sea-level air, $\rho = 1.225$ kg/m³, $g = 9.80665$ m/s².

Find: (a) the drag constant $k$ in kg/m; (b) the terminal velocity $v_t$; (c) the time constant $\tau = v_t/g$; (d) the time to reach 99 % of $v_t$ (use $\operatorname{artanh}(0.99) = 2.6467$); (e) the distance fallen by then. Check (b) with `terminal_velocity_expr().subs(...)` and (e) with `lambdify_position()`. Answer form: five numbers with units, four significant figures.

## 2. (Derivation: solve the equation of motion by hand)

Starting from $m\,\mathrm{d}v/\mathrm{d}t = mg - kv^2$ with $v(0) = 0$:

(a) Separate variables and integrate to obtain $v(t) = v_t\tanh(gt/v_t)$. Name each assumption where you use it (in particular where you need $v < v_t$).
(b) Integrate again to obtain $x(t) = (v_t^2/g)\ln\cosh(gt/v_t)$ with $x(0) = 0$.
(c) Show that for $t \gg v_t/g$, $x(t) \approx v_t t - (v_t^2/g)\ln 2$, and interpret the constant.
(d) Reproduce (a) in SymPy with `sp.dsolve` (extension) or by checking `ode_residual`, and (b) with `sp.integrate`. Report whether SymPy's form is identical to yours and, if not, how you proved them equal.

Answer form: derivation with steps named; one short paragraph for (c); SymPy transcript for (d).

## 3. (Unfamiliar scenario: thrown upward)

The same package (use $m = 1.2$ kg, $k = 0.006125$ kg/m) is launched vertically *upward* at $v_0 = 30$ m/s. Take upward positive; drag now acts downward, and so does gravity.

(a) Write the equation of motion for the upward leg and explain why it is not the same equation as in the lesson with a sign change of $v$.
(b) Solve it (by hand or with SymPy) to get $v(t) = v_t\tan\!\big(\arctan(v_0/v_t) - gt/v_t\big)$ and verify the result with an `ode_residual`-style check.
(c) Find the time to the apex and the apex height, both symbolically and numerically, and compare with the no-drag values $v_0/g$ and $v_0^2/(2g)$.
(d) Why does the lesson's `velocity_expr` not describe the fall *after* the apex if you simply shift the time origin? (It does; say why, in one sentence.)

Answer form: equation, derivation or SymPy transcript, four numbers with units, two sentences.

## 4. (Numerical, uses your code: the error curve)

Run the computational experiment of lesson part 4: `velocity_derivative_check(3.0, 1.2, 9.80665, 0.006125, h)` and `derivative_check(20.0, ...)` for $h = 10^{-1}, \dots, 10^{-8}$.

(a) Tabulate and plot on log-log axes. Fit the slope on each side of the minimum for the velocity check.
(b) Compute the predicted optimum $h_{\text{opt}} = (3\varepsilon|v|/|v'''|)^{1/3}$ at $t_0 = 3$ s using SymPy for $v'''$, and compare with the measured minimum.
(c) Implement the five-point stencil $(-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h))/(12h)$ and repeat (a) for it. Report the new left-hand slope and the new minimum error.

Answer form: a table, two plots, two slopes with uncertainties, one paragraph.

## 5. (Series: when does drag matter?)

(a) From `small_time_series(6)` write down the relative deviation of $v(t)$ from free fall, $1 - v/(gt)$, to leading order in $t$.
(b) Solve for the time $t_1$ at which that deviation reaches 1 % for the lesson's package, then evaluate the *exact* deviation at $t_1$ with `lambdify_velocity` and report the difference.
(c) A P00.15 experiment drops a ball from 2 m and times it to 0.01 s. Using $t_1$ as the yardstick, decide whether free fall is an adequate model for that experiment with a tennis ball ($m = 58$ g, $A = 3.4 \times 10^{-3}$ m², $C_D \approx 0.5$) and with a table-tennis ball ($m = 2.7$ g, $A = 1.26 \times 10^{-3}$ m², $C_D \approx 0.5$).

Answer form: one formula; one time and one percentage; two yes/no decisions each with a number.

## 6. (Judgment: where the model stops)

For each of the following, say which assumption of the lesson's model fails, estimate the size of the effect with one number, and say whether a symbolic closed form still exists:

(a) A 30 km stratospheric balloon drop (ISA density at 30 km is about 0.018 kg/m³).
(b) The package tumbles, with $C_D A$ varying between 0.006 and 0.014 m².
(c) A 20 µm fog droplet falling in air (hint: compute the Reynolds number $\rho v L/\mu$ with $\mu = 1.8 \times 10^{-5}$ Pa·s using the quadratic model's own $v_t$; at $Re < 1$ drag is linear in $v$ and Stokes' law gives $v_t = 2r^2\rho_w g/(9\mu)$). Compare with a 2 mm raindrop.
(d) The fall lasts long enough that $g$ changes by 0.1 % (how far is that?).

Answer form: four short paragraphs, each with a number.

## 7. (Hard, uses your code: linear versus quadratic drag)

A competing model uses **linear drag**, $m\,\mathrm{d}v/\mathrm{d}t = mg - cv$.

(a) Solve it symbolically for $v(0) = 0$ and find its terminal velocity and time constant.
(b) Choose $c$ so that both models have the *same* $v_t$ as the lesson's package. Derive, for each model, the time to reach 90 % of $v_t$ in units of $v_t/g$. Which model reaches terminal speed sooner, and why physically?
(c) Lambdify both $v(t)$ and plot them together with free fall for $0 \le t \le 20$ s. At $t = 2$ s, which is faster and by how much?
(d) Propose one measurement from a P00.15-style drop test that would tell the two models apart, and estimate the precision it needs.

Answer form: two formulas, two dimensionless times, a plot, one number with units, one paragraph.
