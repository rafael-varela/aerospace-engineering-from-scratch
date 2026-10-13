# P00.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and the tolerance you would accept on every numerical answer.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Write it at the top of your work.

## 1. (Warm-up: coefficients in level flight)

A light aircraft of mass $m = 600 + 10s$ kg flies level at $V = 40$ m/s at sea level ($\rho = 1.225$ kg/m³, $g = 9.80665$ m/s²), wing area $S = 12.0$ m². In level flight lift equals weight.

Find: (a) the dynamic pressure $q$ in Pa; (b) $C_L$; (c) $C_D$ and $L/D$ from the default polar $C_{D0} = 0.020$, $k = 0.050$; (d) the drag force in N. Answer form: four numbers with units, each to four significant figures, and a one-line check of (b) with `lift_coefficient`.

## 2. (Derivation: two optima on the polar)

Starting from $C_D = C_{D0} + k C_L^2$:

(a) Derive $C_L^*$ and $(L/D)_{\max}$ by setting $\mathrm{d}(C_L/C_D)/\mathrm{d}C_L = 0$. Show that at the optimum the induced drag equals $C_{D0}$, so $C_D^* = 2C_{D0}$.

(b) The power required for level flight is proportional to $C_D / C_L^{3/2}$. Derive the lift coefficient that minimises power and show it is $\sqrt{3}$ times the lift coefficient of best $L/D$.

(c) State one physical assumption used in (a) that is violated near stall, and say in one sentence what a test for $C_L^*$ would do if the real polar were used instead.

Answer form: derivation with each step named; two formulas; one short paragraph.

## 3. (Floating point by hand)

(a) Convert 0.1 to binary by repeated doubling until the pattern repeats. Write the repeating block.

(b) The spacing between neighbouring `float64` values near a number $x$ is the ulp, $2^{\lfloor \log_2 x \rfloor - 52}$. Compute the ulp near 1531.25 and near 0.3. Check both with `math.ulp`.

(c) `0.5 * 1.225 * 50.0 ** 2` involves two roundings (1.225 is not exact; the product is rounded). Give an upper bound, in ulps and in pascals, on how far the result can be from 1531.25. Then evaluate it in Python and report the actual difference.

(d) Using (c), justify the `rel=1e-12` in `test_dynamic_pressure_analytic_value`. Could it be `rel=1e-15`? Answer form: short calculations and two sentences.

## 4. (Unfamiliar scenario: a wind-tunnel balance in imperial units)

A colleague hands you `coefficient_from_balance(force_lbf, q_psf, area_ft2)`, which should return a dimensionless coefficient from a wind-tunnel balance that logs force in pounds-force and dynamic pressure in pounds-force per square foot. You have not seen the code.

Write, in prose or as Python, **three tests**: one against a known case you can compute by hand, one limiting case, and one invalid-input case. For each, state the tolerance and justify it with a number. Then say which single test would catch a stray factor of 32.174 (the lbm-to-lbf slug conversion) inside the function, and why the others might not.

## 5. (Numerical, uses your code: measure the drift)

Run the computational experiment of lesson part 4 with $m_0 = 50 + s$ kg and a 10 kg total burn, for $n \in \{10, 10^2, 10^3, 10^4, 10^5, 10^6\}$ steps.

Report: the table (`n, dt_s, drift_kg, bound_kg`), the log-log plot, and the slope of $\log(\text{drift})$ against $\log n$ from a straight-line fit over the points where the drift is non-zero. Compare the slope with 1 (the worst-case bound) and 0.5 (random-walk cancellation). Then repeat with the step chosen as exactly 0.125 kg and explain the result. Answer form: table, plot, slope with its uncertainty, two paragraphs.

## 6. (Judgment: when the regression test fails)

You change `K_DEFAULT` from 0.050 to 0.045 because the aircraft got a new wing with higher aspect ratio. `test_polar_table_matches_stored_reference` fails; every other test passes.

(a) Is the code wrong, the test wrong, or neither? What exactly do you do with `REFERENCE_POLAR`, and what do you write in the commit message?

(b) Now suppose instead that you did *not* change anything deliberately, you only "cleaned up" `drag_polar`, and the same test fails with a relative difference of $3 \times 10^{-7}$. List the two most likely causes and the one measurement that distinguishes them.

(c) A reference value for a later lesson will come from `scipy.integrate.solve_ivp` run with `rtol=1e-8`. What is the tightest tolerance you can honestly use in a regression test against it, and why is "tighter is always better" wrong here?

## 7. (Hard, uses your code: does the suite catch planted bugs?)

**Mutation testing** means deliberately breaking the code and checking that at least one test fails. Plant each of the following bugs in a copy of your `main.py`, one at a time, run the suite, and record which tests fail:

1. `drag_polar` returns `cd0 + k * cl` (missing square).
2. `dynamic_pressure` returns `rho * v ** 2` (missing one half).
3. `max_lift_to_drag` returns `(math.sqrt(cd0 / k), 1.0 / math.sqrt(cd0 * k))` (missing the 2).
4. `_check_positive` uses `< 0.0` instead of `<= 0.0`.
5. `is_close` uses `min(abs(a), abs(b))` instead of `max`.

For any bug that *no* test catches, write a new test that does, with its tolerance justified. Answer form: a 5-row table of bug versus failing tests, plus any new test code.
