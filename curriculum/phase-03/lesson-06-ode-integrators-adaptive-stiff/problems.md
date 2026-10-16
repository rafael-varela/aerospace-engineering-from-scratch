# P03.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the stability budget of a stiff model)

Take the lesson's skin-structure model with the skin capacity changed to C_s = 100·s J/K (everything else as in `SkinStructure`).

(a) Estimate the fast eigenvalue by hand from λ_fast ≈ −(k_c/C_s)(1 + C_s/C_b) and the slow one from λ_slow ≈ −k_b/C_b, then compute both with `thermal_jac` and `np.linalg.eigvals` at T = 300 K. Give the stiffness ratio.
(b) Give the largest stable step for explicit Euler and for the Dormand–Prince pair (use `stability_interval`), and the number of steps each needs to cover one hour.
(c) Backward Euler with h = 60 s: what is z = hλ for each mode, and what is |R(z)| for each? One sentence on what each number means for the solution.

Answer form: six numbers with units, two step counts, two |R| values and two sentences.

## 2. (Derivation: the controller, and A- versus L-stability)

(a) Starting from "the local error estimate of a step of size h is E = C h^{q}" (q = 5 for the 5(4) pair), derive the step h_new that would make the estimate exactly equal to the tolerance, and explain in one sentence each why a safety factor below 1 and an upper clamp on the growth are added.
(b) The trapezoidal rule is y_{n+1} = y_n + (h/2)[f(t_n, y_n) + f(t_{n+1}, y_{n+1})]. Derive its stability function R(z), show that |R(z)| < 1 for every z with Re z < 0 (so it is A-stable), and find lim_{z→−∞} R(z). Compare with backward Euler and explain why the limit matters for a mode with hλ = −300.
(c) For the linear two-node model (ε = 0) write the 2×2 matrix A, its trace and determinant, and show that the product of the eigenvalues is k_c k_b/(C_s C_b). Use trace and determinant to recover λ_fast and λ_slow for the lesson's parameters without an eigen-solver, and compare with the computed −5.0223 and −1.2414 × 10⁻⁴ s⁻¹ (which include the radiation term; estimate its contribution).

Answer form: three short derivations with the formulas boxed.

## 3. (Unfamiliar scenario: an actuator in the loop of a flight simulation)

A longitudinal flight-dynamics model has a phugoid mode with period 60 s and a short-period mode with period 3 s (both lightly damped). An elevator actuator is added as a first-order lag with time constant τ = s milliseconds. The simulation must run 10 minutes of flight.

(a) Write the eigenvalue the actuator contributes and the approximate eigenvalues of the two aircraft modes (imaginary parts from the periods; ignore damping). What is the stiffness ratio of the combined system?
(b) Estimate the step the Dormand–Prince pair must keep for stability and the number of steps and right-hand-side evaluations for 10 minutes. Compare with the step that *accuracy* on the short-period mode would need at rtol = 10⁻⁶ (take the rule of thumb "about 40 steps per period for a fifth-order method at 10⁻⁶").
(c) An engineer proposes to "fix" the cost by replacing the actuator lag with its steady-state value (δ_e follows the command instantly). Is the resulting system still stiff? What is lost?
(d) Which `solve_ivp` method would you recommend if the actuator must stay, and what would you supply to it besides the right-hand side?

Answer form: three eigenvalues, one ratio, two step counts, two short paragraphs.

## 4. (Numerical, uses your code: your own entry angle)

Set the entry flight-path angle to γ_E = −(4 + s/10) degrees (everything else as in `entry_state`).

(a) From (7), predict the Allen–Eggers peak deceleration in g and the altitude h* of the peak.
(b) Run `solve_rk45` on the full model (gravity on) with atol = 10⁻³ at rtol = 10⁻⁶ and 10⁻⁸, stopping below 20 km. Report accepted steps, rejections, evaluations, and the peak deceleration and altitude read from the stored points.
(c) Run once more at rtol = 10⁻¹⁰ with `h_max = 0.25` s to get a finely resolved peak. How far are the (b) peaks from it, and why does the rtol = 10⁻⁸ run read a *lower* peak than the 10⁻⁶ run?
(d) Compute the step-count growth exponent log(N₈/N₆)/log 100 and compare with the 1/5 of the theory. Which part of the integration (start-up, peak, tail) breaks the assumption behind the theory?

Answer form: a table with two rows (rtol, steps, rejections, nfev, peak g, h_peak), four numbers from (a) and (c), a short paragraph.

## 5. (Numerical, uses your code: a second-order implicit method)

Write `trapezoidal_step(f, t, y, h, jac)` with the same Newton structure as `implicit_euler_step` (start from z = y; G(z) = z − y − (h/2)[f(t, y) + f(t + h, z)]; G′ = I − (h/2)J).

(a) Verify on y′ = λy that one step gives y·(1 + hλ/2)/(1 − hλ/2).
(b) Run the thermal model from T_s(0) = 400 K, T_b(0) = 300 K (the skin starts 100 K hot) for one hour with h = 60, 30 and 15 s, with both trapezoidal and backward Euler, against a Radau reference (rtol 10⁻¹²). Print the first four skin temperatures of each run and the final errors.
(c) Explain what you see in the trapezoidal skin temperatures using the limit from Problem 2(b). Why is the observed order of the trapezoidal rule across these three runs not 2?
(d) Repeat the trapezoidal convergence study at h = 1, 0.5 and 0.25 s and report the observed order. What is the condition for the asymptotic order to appear?

Answer form: code, two tables, two short paragraphs with the number |R(−301)| in them.

## 6. (Judgment: choosing the integrator for a coupled simulation)

You are asked to simulate the entry trajectory of part 2 *together with* the skin-structure thermal model, with the heat input q set by the trajectory (q ∝ ρ^{1/2} v³, the usual convective-heating scaling), for 200 s of flight, in a Monte Carlo study of 10 000 dispersed cases that must finish overnight on one laptop core.

(a) Which part of the coupled system sets the stable step of an explicit method, and what is that step?
(b) Estimate the right-hand-side evaluations per case for (i) adaptive RK45, (ii) backward Euler with h = 1 s, (iii) `solve_ivp` with `Radau` at rtol 10⁻⁶, using the step counts from this lesson as calibration. Which fit the overnight budget if one evaluation costs 20 μs?
(c) The trajectory needs rtol = 10⁻⁸ for the dispersion statistics to be meaningful, but the thermal model's inputs are uncertain to 20 %. How do you reconcile the two tolerances inside one solver call?
(d) A colleague suggests integrating the trajectory with RK45 and the thermal model separately with backward Euler, exchanging q and the state once per second. Name one thing this gains and one thing it risks, with a sentence on how you would test whether the risk matters.

Answer form: four short paragraphs, each with at least one number.
