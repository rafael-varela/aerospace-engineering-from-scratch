# P02.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical (and with what step or tolerance).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a burn by hand)

A motor burns a 300 kg rocket down to 120 kg at ṁ = 5 + s/10 kg/s with exhaust velocity 2400 m/s, vertically, in vacuum, with g = 9.81 m/s².

(a) Find the burn time, the Tsiolkovsky Δv, the gravity loss and the burnout velocity by hand from (1)–(3). (b) Find the burnout altitude from (4). (c) Find the vacuum apogee. (d) Confirm (a) and (b) with `vacuum_velocity`, `vacuum_altitude` and with `simulate_rocket(..., drag=False)` using 400 RK4 steps; report the RK4 error.

Answer form: five numbers with units, one RK4 error.

## 2. (Derivation: the altitude integral and the drag coast)

(a) Starting from v(t) = v_e ln(m₀/m(t)) − gt with m(t) = m₀ − ṁt, derive h(t) = v_e t − (v_e m(t)/ṁ) ln(m₀/m(t)) − ½gt² by substitution and integration by parts. Show every step, including the check h(0) = 0.
(b) For the upward coast dv/dt = −g − kv² (constant k), derive the height gained to the apex, Δh = ln(1 + kv_b²/g)/(2k), by writing dh = v dv/(dv/dt) and integrating in v. State where the sign of the drag term would change and why the formula stops there.
(c) Expand Δh for small k to two terms and interpret the correction term physically. Verify the expansion numerically with `quadratic_drag_coast` at k = 10⁻⁶ and 10⁻⁵ m⁻¹ for v_b = 1000 m/s.

Answer form: two derivations, one expansion with a numerical check.

## 3. (Unfamiliar scenario: an avionics box warming up)

An avionics box of heat capacity C = 2000 + 50s J/K is switched on at −10 °C after a cold soak. Its electronics dissipate P = 30 W and it loses heat to the 20 °C cabin through a conductance G = 1.5 W/K, so C dT/dt = P − G(T − T_a).

(a) Classify the equation (order, linear or not, autonomous or not, separable or not) and write it in the standard first-order linear form. (b) Solve it with an integrating factor and identify the time constant and the steady-state temperature. (c) How long until the box reaches its 15 °C operating minimum? (d) Check (c) by integrating the equation with `integrate_fixed_step` (RK4, 200 steps) and by `equilibria` on the autonomous form; state what the phase line says about T → T_∞.

Answer form: a classification, a solution with two named constants, a time in minutes, two checks.

## 4. (Numerical, uses your code: your rocket's drag)

Take the sounding rocket of the lesson but with C_D = 0.3 + s/100.

(a) Run `simulate_rocket` with RK4 to 150 s and report the burnout velocity and altitude, the apogee and the time of apogee. (b) Repeat with 300, 600, 1200 and 3000 steps and show the apogee has converged to 1 m. (c) Compute the frozen-density coast bracket (`quadratic_drag_coast` with ρ at the burnout altitude) and the vacuum coast, and show the full-model height gain lies between them. (d) Compare the burnout velocity with `solve_ivp` at rtol = 10⁻¹⁰.

Answer form: a table of four states, a convergence table, two bracket values, one oracle difference.

## 5. (Judgment: the way down)

The empty rocket (100 kg, C_D = 0.5, A = 0.0707 m²) falls back from a 45 km apogee.

(a) Compute its terminal velocity at 45, 30, 20, 10 and 0 km with `terminal_velocity` and the exponential atmosphere. (b) Integrate the descent with `simulate_rocket` (run the baseline to 400 s) and compare the actual speed at 18 km and at the ground with the local terminal velocities. Is the descent quasi-steady (does v track the local terminal velocity)? (c) Explain the answer to (b) with the time constant of the autonomous equation near its equilibrium, 1/|f′(v*)| = 1/(2k|v*|), compared with the time the rocket spends in each density layer. (d) Name the assumption of the descent model that fails first for a real rocket body and one measurable sign that it has.

Answer form: a table of five terminal velocities, two comparisons, a time-constant argument, a judgment.

## 6. (Hard: the global error of Euler's method)

(a) Show that for dy/dt = f(t) (no y on the right) Euler's method is the left Riemann sum of ∫f dt and that its global error on [a, b] is ½h[f(b) − f(a)] + O(h²). Check it against the Euler error constant the lesson reports for the vacuum burn.
(b) For dy/dt = f(t, y) with |∂f/∂y| ≤ L and |y″| ≤ M, show by induction on the recursion e_{n+1} ≤ (1 + hL)e_n + ½h²M that the global error satisfies e_n ≤ (hM/2L)(e^{L(t_n − t₀)} − 1).
(c) Apply the bound to dv/dt = −g − kv|v| with k = ρ₀C_D A/(2m) for the empty rocket at sea level, from rest, over 60 s, and compare with the measured Euler error at h = 1, 0.5, 0.25, 0.125 and 0.0625 s (use `solve_ivp` at rtol = 10⁻¹² as the reference). By what factor is the bound pessimistic, and why?

Answer form: two derivations, one table with a ratio.
