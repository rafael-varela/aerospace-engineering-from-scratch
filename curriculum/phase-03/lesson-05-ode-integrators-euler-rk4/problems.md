# P03.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are bounds and which are observed.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one step by hand, then many by code)

Take y′ = −y + t, y(0) = 1, whose exact solution is y = t − 1 + 2e^{−t}, and the step h = 0.1 + s/100.

(a) Compute one step by hand with forward Euler, Heun and RK4, writing out every kᵢ. Compare each with y(h) and report the three errors.
(b) With the code, integrate to t = 2 with h, h/2 and h/4 for the three methods. Report the errors at t = 2, the two consecutive ratios for each method, and the observed orders.
(c) Why are the ratios not exactly 2, 4 and 16?

Answer form: the three hand steps with their k's, a 3 × 3 error table, six ratios, three orders, two sentences for (c).

## 2. (Derivation: local and global error, the stability function, RK4 as Simpson)

(a) Show that one forward Euler step from the exact solution has error h²y″(ξ)/2, and that after n = T/h steps on a problem with Lipschitz constant L the global error is at most (e^{LT} − 1)/L × (h/2) max|y″|. Hint: let eₙ be the error after n steps and show eₙ₊₁ ≤ (1 + hL)eₙ + h²M/2.
(b) Show that Heun's method applied to y′ = λy multiplies y by 1 + z + z²/2 with z = hλ, and that RK4 multiplies it by the series through z⁴/24. Deduce the real-axis stability limits of Heun (exactly 2) and RK4 (numerically 2.7853), and show that Euler is unstable for every h on y′ = iωy.
(c) Show that on y′ = g(t) the RK4 step is Simpson's rule for ∫g over [t, t + h], hence exact when g is a cubic.
(d) Explain in two sentences why a one-step method with local error O(h^{p+1}) has global error O(hᵖ), and why the e^{LT} factor in (a) is unavoidable in general but pessimistic for a dissipative system.

Answer form: four derivations.

## 3. (Unfamiliar scenario: apogee of a sounding rocket in coast)

After burnout a 20 kg sounding rocket coasts vertically with velocity v₀ = (200 + 5 s) m/s through air of density 1.225 kg/m³ with drag area C_D A = 0.02 m² (SYNTHETIC: constant density and drag coefficient so that a closed form exists). The equations are ż = v, v̇ = −g − (ρC_DA/2m) v|v| with g = 9.80665 m/s². For the upward leg the exact solution is v(t) = v_t tan(arctan(v₀/v_t) − gt/v_t) with terminal speed v_t = √(2mg/(ρC_DA)), the apogee time is t_a = (v_t/g) arctan(v₀/v_t), and the apogee height is z_a = (v_t²/2g) ln(1 + v₀²/v_t²).

(a) Compute v_t, t_a and z_a. Compare z_a with the no-drag value v₀²/2g.
(b) Integrate with Euler, Heun and RK4 at h = 0.5, 0.25 and 0.125 s to t = 14 s (before apogee) and report the altitude error at t = 14 s and the observed order of each method.
(c) Integrate past apogee (to 16 s), locate the apogee from the sign change of v (linear interpolation of v in the bracketing step, then z from a parabola through the bracketing state), and report the apogee time and height errors for the three methods at the three steps. Which converge at the method's order?
(d) A flight computer must predict apogee in real time from this model with a 0.1 m error. Which method and step would you use, and how many evaluations of the right-hand side does a 15 s prediction cost?

Answer form: three numbers for (a) with the comparison, a 3 × 3 table and three orders for (b), two 3 × 3 tables for (c) with a one-line verdict, a method, a step and a count for (d).

## 4. (Uses the code: a different orbit)

Take a = (8000 + 200 s) km and e = 0.1 + s/100, μ = 3.986004418 × 10¹⁴ m³/s².

(a) Compute the period, the perigee altitude and the perigee and apogee speeds from `kepler_state`.
(b) Tabulate the one-orbit position error of Heun, Kutta3 and RK4 at 200, 400 and 800 steps per orbit and give the observed orders.
(c) From the 800-step RK4 error and formula (2) predict the steps per orbit for a 1 m error after one orbit; check with a run.
(d) Measure the RK4 error after 1, 2 and 4 orbits at 800 steps per orbit. Assuming the growth law you observe, predict the steps per orbit that keep the error below 1 m after one day (86 400 s), check it with a run, and give the evaluation count for the day.
(e) Report the relative energy drift of RK4 and Heun after five orbits at 400 steps per orbit, with signs.

Answer form: four numbers for (a); a 3 × 3 table and three orders for (b); a prediction and a check for (c); the growth law, a prediction, a check and a count for (d); two signed numbers for (e).

## 5. (Uses the code: stiffness and the stability limit)

Take y′ = −λ(y − cos t), y(0) = 0, with λ = 100 + 10 s. The exact solution is y = (λ² cos t + λ sin t)/(1 + λ²) − λ²e^{−λt}/(1 + λ²): after a transient of a few 1/λ it is cos t to within 1/λ.

(a) Compute h_max for Euler, Heun and RK4 and verify with runs at 0.95 h_max and 1.05 h_max to t = 10.
(b) How many evaluations of f does each method need to reach t = 10 at 0.95 h_max? What is RK4's error at t = 10 there, and what step would RK4 need for a 10⁻⁶ error if stability were not an issue (use (2) from two runs at h = 0.01 and 0.005)?
(c) Call `scipy.integrate.solve_ivp` with `RK45`, `Radau` and `LSODA` at rtol = 10⁻⁶ and compare `nfev`. Which methods are hostage to stability?
(d) Now take the undamped oscillator y″ = −y. Show from the stability function that Euler's amplitude grows by (1 + h²)^{n/2} over n steps, evaluate the growth per period for h = 0.1 and 0.01, and compute |R(ih)| for RK4 at the same steps. What does this predict for a numerical orbit?

Answer form: three h_max with six run results for (a); three counts, an error and a step for (b); three nfev with a sentence for (c); two growth factors, two RK4 moduli and a two-sentence prediction for (d).

## 6. (Judgment: a method and a step for a one-day propagation)

A conjunction-screening tool must propagate 2000 catalogue objects, all in orbits like the lesson's (a ≈ 7500 km, e ≤ 0.1), for one day with a position error below 10 m per object relative to the two-body truth, on a machine that evaluates the two-body right-hand side 2 × 10⁶ times per second in NumPy.

(a) From the lesson's measured errors (4.7 mm at 800 RK4 steps per orbit, 4.8 km for Heun; linear growth in time), pick a method and a step for 10 m after 13.4 orbits and give the evaluation count per object and the total run time.
(b) Would Euler or Heun ever be the right choice here? Give the evaluation counts.
(c) The catalogue also contains objects with e = 0.7. What goes wrong with the fixed step, and what do you tell the tool's author to do about it (two sentences)?
(d) The tool's author proposes to validate the integrator by comparing with `solve_ivp` at default tolerances (rtol = 10⁻³). Why is that the wrong oracle, and what would you compare with instead?

Answer form: four short paragraphs; (a) and (b) must contain counts.

## 7. (Proof: the order conditions of two-stage methods and the stability function)

Consider the general explicit two-stage method k₁ = f(t, y), k₂ = f(t + c₂h, y + c₂hk₁), yₙ₊₁ = y + h(b₁k₁ + b₂k₂) for the autonomous scalar problem y′ = f(y).

(a) Expand the exact solution to h³: y(t + h) = y + hf + h²f′f/2 + h³(f″f² + f′²f)/6 + O(h⁴), where f′ = df/dy.
(b) Expand k₂ to h² and show the method's step is y + h(b₁ + b₂)f + h²b₂c₂f′f + h³b₂c₂²f″f²/2 + O(h⁴).
(c) Deduce the order-2 conditions (1), and show that no choice of b₁, b₂, c₂ gives order 3 (compare the two h³ terms).
(d) Show that for any explicit s-stage method applied to y′ = λy the stability function R(z) is a polynomial of degree ≤ s, and that if the method has order p then R(z) = Σ_{k≤p}z^k/k! + O(z^{p+1}). Conclude that for s = p ≤ 4 the stability function is exactly the truncated exponential, so Euler, both RK2 methods, Kutta3 and RK4 have the functions used in the lesson.

Answer form: four derivations.
