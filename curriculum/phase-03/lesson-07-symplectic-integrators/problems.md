# P03.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: Verlet on the oscillator at your own step)

Take q'' = −ω²q with hω = s/20.

(a) Write the 2×2 Verlet map M of equation (3) with numbers, and check det M = 1.
(b) Give the coefficient (1 − h²ω²/4) of the modified energy (4) and the relative amplitude of the true energy's oscillation.
(c) Give the numerical frequency ω_h/ω = (2/(hω)) asin(hω/2) and its relative error; compare with the leading term h²ω²/24. How many steps per period is your hω, and how many is the stability limit?

Answer form: a matrix, four numbers, one sentence.

## 2. (Derivation: three properties of kick-drift-kick)

(a) Show that a kick p ← p + h a(q) and a drift q ← q + h p each have a symplectic Jacobian (MᵀJM = J) when ∂a/∂q is symmetric, and hence that symplectic Euler, Verlet and Yoshida's composition are symplectic for any force derived from a potential.
(b) Show that for a central force a(q) = f(|q|) q, a kick and a drift each leave L = q × p unchanged, so these methods conserve angular momentum exactly at any step.
(c) Show that Verlet is symmetric: applying the step with −h to (q_{n+1}, p_{n+1}) returns (q_n, p_n) exactly. Then show that symplectic Euler is not symmetric, and that N symplectic Euler steps equal N drift–kick–drift Verlet steps conjugated by a half-step drift, SE_h^N = D_{h/2} ∘ Ṽ_h^N ∘ D_{−h/2} (which is why its position error over a whole period scales as h²).

Answer form: three short proofs.

## 3. (Unfamiliar scenario: a Molniya orbit with a fixed step)

A communications satellite is in a Molniya-type orbit: a = 26 600 km, e = 0.74 (perigee 6 916 km, apogee 46 284 km, period 12.0 h).

(a) Compute the perigee and apogee speeds and the angular rates ω_p = v_p/r_p and ω_a = v_a/r_a. Show that ω_p/ω_a = ((1 + e)/(1 − e))².
(b) To keep hω_p ≤ 0.05 at perigee with a fixed step, how many steps per orbit are needed? How many would keep hω_a ≤ 0.05 at apogee? What fraction of the first count is "wasted" on the slow part of the orbit?
(c) Run Verlet and RK4 for 5 orbits at 2 000 and 8 000 steps per orbit and report the maximum |ΔE/E|. Which method would you use for 1 000 orbits of this satellite at a fixed step, and why does the answer differ from the e = 0.1 case of the lesson?
(d) Propose, in two sentences, a change of independent variable that would make the step small at perigee and large at apogee without an error estimate (extension 1), and say what it must preserve.

Answer form: six numbers, one table, two short paragraphs.

## 4. (Numerical, uses your code: your own eccentricity)

Keep the perigee at 6 300 km and set e = 0.1 + s/100 (so a = 6 300 km/(1 − e)). Propagate 50 orbits at 200 steps per orbit with Verlet, RK4 and Yoshida.

(a) Give a, the period, v_p, ω_p and hω_p.
(b) Report for each method the maximum |ΔE/E|, the final ΔE/E, the fitted slope per orbit, the maximum |ΔL/L| and the final position error against `kepler_exact`.
(c) Predict the ratio of Verlet's energy swing to the e = 0.1 value (1.06 × 10⁻⁴) from (hω_p)² alone, and compare with what you measure. Why is the measured ratio larger?
(d) Is RK4's final energy error or Verlet's larger? Which would you rather have after 5 000 orbits, and what number from your table supports the choice?

Answer form: five numbers, a table of three rows, a short paragraph.

## 5. (Numerical, uses your code: matched cost and the crossover)

Verlet with force reuse costs one evaluation per step; RK4 costs four. Compare Verlet at 800 steps per orbit with RK4 at 200 (equal cost) on the e = 0.1 orbit.

(a) Run Verlet for 20 orbits and record its energy swing; run RK4 for 100 orbits and fit its drift slope. Compute the crossover orbit count N* = swing/|slope| at which RK4's secular error first exceeds Verlet's bound.
(b) Run Yoshida at 200 steps per orbit (600 evaluations) for 20 orbits and record its swing. Where does it sit?
(c) RK4's position error after 60 orbits is 1.0 km and Yoshida's 2.5 km at the same step. Reconcile this with (b): what is each number measuring?
(d) For a 10⁴-orbit propagation at 800 evaluations per orbit, which of the three would you run, and what is the predicted |ΔE/E| at the end for each?

Answer form: four numbers, one short paragraph, one recommendation with three numbers.

## 6. (Judgment: when symplectic is the wrong tool)

For each situation say whether you would use a symplectic fixed-step method, an adaptive Runge–Kutta method, or either, and give the deciding property in one sentence: (a) a lunar transfer with two finite burns and a flyby, 5 days; (b) a 30-year propagation of 10 000 debris fragments above 1 000 km altitude; (c) the re-entry of `P03.06`; (d) a 10⁶-year integration of the outer planets; (e) a satellite at 400 km with drag for one year; (f) a molecular-dynamics simulation of 10⁵ atoms for 10⁶ steps.

Answer form: six labelled lines.
