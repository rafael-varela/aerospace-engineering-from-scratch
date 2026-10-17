# P03.02  Root Finding

> Where is the satellite at 09:30? Kepler wrote the equation in 1609 and nobody has solved it in closed form since. Four iterations do, and the difference between them is how many correct digits each step buys.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the convergence of bisection (rate ½), fixed-point iteration (rate |g′(x*)|), Newton's method (order 2 at a simple root, rate ½ at a double root) and the secant method (order (1 + √5)/2) from Taylor expansions, and state the conditions each needs. 2. Implement and test all four from scratch with their full iterate histories, measure the observed convergence order from those histories, and match it to theory. 3. Solve Kepler's equation M = E − e sin E for any mean anomaly and eccentricity, convert to true anomaly and position–velocity, and verify the result against vis-viva, angular momentum and `scipy.optimize.brentq`. 4. Decide which method a flight computer should use given its tolerance, its cost budget, whether a derivative is available and how far the starting guess may be from the root. |
| **Prerequisites** | `P03.01` (rounding floor, observed order), `P01.05` (Taylor expansions with remainder; the small-angle series), `P02.05` (the Clohessy–Wiltshire orbit constants n, μ, a), `P01.02` (the Newton step as the tangent-line zero). Skill check: you can write f(x + h) = f(x) + h f′(x) + h² f″(ξ)/2 and compute one Newton step by hand. |
| **Engineering question** | A satellite is in a geostationary transfer orbit (perigee 250 km, apogee 35 786 km, period 10.53 h). Where is it, in position and velocity, 30 and 60 minutes after perigee, and how long after perigee does it cross 22 500 km radius? Each answer requires the eccentric anomaly E from Kepler's equation, which has no closed form: which iteration solves it to 10⁻¹⁴ rad in the fewest function evaluations, and which one is guaranteed never to fail? |
| **Assumptions and limits** | Two-body Keplerian motion (point-mass Earth, no J₂, drag or third bodies), elliptical orbits 0 ≤ e < 1 (the hyperbolic case is Problem 3), angles in radians with M reduced to one revolution; f continuous for bisection, C² with f′(x*) ≠ 0 for the Newton and secant orders; the stopping test |x_{k+1} − x_k| ≤ tol, which for a fast method is a safe proxy for the error but for a slow one (fixed point at rate 0.9) understates it by a factor 1/(1 − rate). The convergence-order measurement works only on the iterates above the rounding floor of `P03.01`. |
| **Analytical work** | The error recurrences e_{k+1} ≤ e_k/2 (bisection), e_{k+1} ≈ g′(x*) e_k (fixed point), e_{k+1} ≈ (f″/2f′) e_k² (Newton), e_{k+1} ≈ (f″/2f′) e_k e_{k−1} and q² = q + 1 (secant); the step count ⌈log₂((b − a)/tol)⌉ for bisection; the contraction constant e of the Kepler fixed-point map; the derivative 1 − e cos E ∈ [1 − e, 1 + e]; the conversions E → ν → (r, v) and the invariants they must satisfy. |
| **Computational deliverable** | `starter/main.py`: `bisection`, `bisection_steps_needed`, `newton`, `secant`, `fixed_point`, `errors_from_iterates`, `observed_convergence_order`, `fitted_convergence_order`, `observed_linear_rate`, `kepler_residual`, `kepler_derivative`, `solve_kepler`, `true_anomaly`, `mean_motion`, `orbit_state`, `specific_energy`, `rayleigh_pitot_ratio`, `mach_from_pitot_ratio`. Tests prove: bisection takes exactly ⌈log₂((b − a)/tol)⌉ midpoints and every midpoint obeys its bound; Newton reaches √2 in ≤ 6 steps with per-step order 2.0 ± 0.1 and equals `scipy.optimize.newton`; Newton raises on a zero derivative, on divergence (atan from 1.5) and on the 0 → 1 → 0 cycle of x³ − 2x + 2; at a double root its rate is exactly ½; the secant method equals `brentq` and raises on a flat slope; fixed-point iteration on cos x converges at rate |sin x*| = 0.674 and diverges for g′ = 2; on Kepler's equation the fitted orders are 2.0 ± 0.2, 1.618 ± 0.2 and 1.00 ± 0.05 for Newton, secant and fixed point, and Newton < secant < bisection < fixed point in steps at e = 0.9; all four methods agree with `brentq` to 10⁻¹² for five eccentricities up to 0.99 and nine mean anomalies, with residuals below 10⁻¹³; revolutions are carried; the orbit state satisfies vis-viva, |r × v| = √(μa(1 − e²)) and ε = −μ/2a to 10⁻¹² at twelve random times and returns to perigee after one period to a millimetre; the Rayleigh pitot inversion round-trips to 10⁻⁹ in Mach number. |
| **Practical activity** | **computational experiment**: for one orbit's eccentricity and mean anomaly, run the four methods from the same start, tabulate error against iteration and function evaluations, measure the order of each, then sweep e toward 1 and M toward 0 to find where Newton's step count climbs and where the fixed-point map becomes slower than bisection. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `root_finding.png`, the convergence table (method, steps, evaluations, fitted order, theory), the e-sweep table, the hand derivations of the four error recurrences, and a half-page note recommending a Kepler solver (method, starting guess, stopping rule, iteration cap) for a flight computer and saying what it does when e → 1. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the Newton and secant recurrences are derived and the secant order obtained from q² = q + 1; all public tests pass; the measured orders are reported with the number of usable iterates; the note states a method, a cap and a fallback. |
| **References** | `@ocw-18-330` "Lecture notes 9–10: Nonlinear Equations" (required-free): bisection, Newton, secant and their convergence orders. `@ocw-18-335j` Week 1 "Lecture 1: Course Overview, Newton's Method for Root-Finding" (handout "Square Roots via Newton's Method") (required-free): quadratic convergence seen digit by digit. `@scipy-userguide` "Optimization (scipy.optimize)" > "Root finding" > "Scalar functions" (required-free): `root_scalar`, `brentq`, `newton` and when to bracket. |
| **Extensions** | ★ Implement Halley's method (third order, uses f″) for Kepler's equation and count evaluations against Newton at the same tolerance. ★ Implement Brent's method (bisection + secant + inverse quadratic interpolation) and compare its step count with `scipy.optimize.brentq` on the Kepler grid. ★ Build a starting-guess table E₀(M, e) that keeps Newton under four steps for every e ≤ 0.99. |
| **Time** | unmeasured |

## 1. The engineering problem

A geostationary transfer orbit (GTO) has perigee 250 km above the Earth's surface and apogee at the geostationary radius 42 164 km: semi-major axis a = 24 396 km, eccentricity e = 0.7283, period T = 10.53 h. Mission control needs to know where the satellite is at any time t after perigee: to point an antenna, to plan the apogee burn, to say when the spacecraft crosses the 22 500 km radius where a ground station loses it. Kepler's laws give the answer in two explicit steps and one implicit one. The **mean anomaly** M = n t, with n = √(μ/a³) the mean motion, grows uniformly in time. The position on the ellipse is parameterised by the **eccentric anomaly** E through r = a(cos E − e, √(1 − e²) sin E). And the two are linked by **Kepler's equation**

**M = E − e sin E**,   (1)

which is transcendental: no finite combination of elementary functions gives E from M. Since 1609, every orbit computation has solved (1) numerically, and the question of this lesson is how.

The question matters beyond orbits. Every "find the x where f(x) = 0" in engineering (the trim angle of attack where pitching moment vanishes, the Mach number behind a pitot reading, the altitude where density reaches a value, the time when a range equals a threshold, the eigenvalue where a determinant vanishes) is the same problem, and the same four iterations answer it. They differ in three ways that an engineer must trade: **how fast they converge** (how many correct digits each step adds), **what they cost** (function evaluations, derivative evaluations), and **when they fail** (no bracket, bad starting guess, zero derivative). A flight computer that must solve (1) at 10 Hz for e = 0.73 in a fixed time budget wants the fastest method that cannot fail; this lesson gives it the tools to make that choice with numbers.

Success criteria, written first: every method solves (1) to 10⁻¹⁴ rad and agrees with `scipy.optimize.brentq`; the measured order of each matches its theory (2, 1.618, 1, 1); the resulting position and velocity satisfy the orbit's conservation laws to 10⁻¹²; and the lesson ends with a recommendation that names a method, a start, a stopping rule and a cap.

## 2. Understand and derive

Throughout, x* is a root of f, e_k = |x_k − x*| is the error of the k-th iterate, and f is as smooth as each method needs. 18.330 notes 9–10 cover the same ground; the derivations below are the ones the tests check.

### 2.1 Bisection: guaranteed and slow

If f is continuous and f(a) f(b) < 0, the intermediate value theorem puts a root in [a, b]. Evaluate the midpoint m, keep the half where the sign changes, repeat. After k midpoints the bracket has width (b − a)/2ᵏ and the k-th midpoint is within (b − a)/2ᵏ of a root, so

**e_k ≤ (b − a)/2ᵏ**, and k = ⌈log₂((b − a)/tol)⌉ midpoints guarantee e_k ≤ tol.   (2)

Each step costs one evaluation and gains exactly one binary digit: 48 steps for a bracket of width 1.46 rad to 10⁻¹⁴. This is **linear convergence with rate ½**: the error *bound* halves every step. The error of the midpoint itself jumps around inside the bound (a midpoint can land very near the root and the next one is far again), which is why the per-step order estimates of 2.3 are meaningless for bisection and the tests check the bound instead. Bisection needs no derivative, no smoothness beyond continuity and no good starting guess; it cannot diverge. What it cannot do is go fast, find a root of even multiplicity (no sign change), or work on vector problems.

### 2.2 Fixed-point iteration: as fast as the slope allows

Rewrite f(x) = 0 as x = g(x) and iterate x_{k+1} = g(x_k). If x* = g(x*) and g is C¹, the mean value theorem gives x_{k+1} − x* = g(x_k) − g(x*) = g′(ξ)(x_k − x*), so

**e_{k+1} ≈ |g′(x*)| e_k**:   (3)

linear convergence with rate |g′(x*)| if that is below 1, divergence if above. For Kepler's equation the natural map is g(E) = M + e sin E with g′(E) = e cos E, so |g′| ≤ e < 1 everywhere: the map is a **contraction** with constant e and converges from any start (Problem 7 proves it). The rate at the root is e|cos E*|: for e = 0.73 and M = 60° it is 0.1399 and the iteration takes 17 steps to 10⁻¹⁴ (observed rate 0.1399); for e = 0.9 and M = 10° it is 0.886, and 63 steps are needed, more than bisection. The number of iterations scales as log(tol)/log(rate), and a rate of 0.99 means 3200 steps: fixed-point iteration is the method to use when g′ is small and to avoid otherwise.

### 2.3 Newton's method: digits that double

Expand f about x_k: 0 = f(x*) = f(x_k) + (x* − x_k) f′(x_k) + (x* − x_k)² f″(ξ)/2. Dropping the quadratic term and solving gives the **Newton step** x_{k+1} = x_k − f(x_k)/f′(x_k), the zero of the tangent line (`P01.02`). Keeping it gives the error:

**e_{k+1} = |f″(ξ)/(2f′(x_k))| e_k² → (|f″(x*)|/2|f′(x*)|) e_k²**.   (4)

The error is *squared* each step: **quadratic convergence**, order 2. Eight correct digits become sixteen in one step, and the method reaches the rounding floor from e₀ = 0.08 in four steps for Kepler at e = 0.73 (errors 8.4 × 10⁻², 2.4 × 10⁻³, 1.8 × 10⁻⁶, 9.9 × 10⁻¹³, 2.2 × 10⁻¹⁶; per-step order estimates 2.02, 2.00). The constant f″/2f′ sets how early the squaring starts: for (1), f″ = e sin E and f′ = 1 − e cos E, and near perigee with e → 1 the constant is large and the first steps can overshoot wildly (for e = 0.99, M = 5° the iterates in degrees run 9.9, 205, 87, 60, 48, 45.5, 45.36: three wasted steps before the squaring takes hold, nine steps in all).

Three failure modes, each in the tests. **Zero derivative**: the tangent is horizontal and the step is infinite (x² − 1 from x₀ = 0). **Divergence**: for f = atan x from x₀ = 1.5 each step overshoots further (|x₀| > 1.39 is enough) and the iterates run to 10¹⁵⁰ in a few dozen steps. **Cycles**: for x³ − 2x + 2 from x₀ = 0 the iteration visits 0, 1, 0, 1, … for ever. **Multiple roots** are a fourth: at a double root f′(x*) = 0, (4) breaks down, and the recurrence becomes e_{k+1} = e_k/2 exactly for f = (x − 1)²: Newton degrades to linear convergence at rate ½, 34 steps instead of 6 (Problem 2b derives this). A production Newton therefore caps its iterations, checks |x_{k+1} − x_k| is shrinking, and falls back to a bracketing method when it is not. For Kepler's equation none of these can happen: f′ = 1 − e cos E ≥ 1 − e > 0, f is monotone with one root per revolution, and the starting guess E₀ = M + e sin M (one fixed-point step from M) keeps Newton under ten steps even at e = 0.999.

### 2.4 Secant method: Newton without the derivative

Replace f′(x_k) by the slope through the last two points: x_{k+1} = x_k − f(x_k)(x_k − x_{k−1})/(f(x_k) − f(x_{k−1})). One new evaluation per step and no derivative. The error recurrence (Problem 2d) is

**e_{k+1} ≈ (|f″|/2|f′|) e_k e_{k−1}**.   (5)

Seeking e_{k+1} = C e_k^q gives e_k e_{k−1} ∝ e_k · e_k^{1/q}, so q = 1 + 1/q, **q² = q + 1, q = (1 + √5)/2 = 1.618**: superlinear, between the fixed point's 1 and Newton's 2. Per *evaluation* the secant method is often the fastest of the four: two secant steps cost two evaluations and multiply the digits by 1.618² = 2.6, one Newton step costs two evaluations (f and f′) and doubles them. Its per-step order estimates oscillate (4.4, 1.2, 1.8 on the Kepler example) because (5) involves two previous errors, not one; the least-squares fit of log e_{k+1} against log e_k is the steadier measurement (1.53 to 1.71 on the Kepler cases in the tests, 1.63 on x³ − 2). The secant method fails when f(x_k) = f(x_{k−1}) (flat slope) and, like Newton, can diverge from a poor start.

### 2.5 Kepler's equation and what to do with E

Reduce M to [0, 2π) (remember the whole revolutions), solve (1), add them back. The derivative 1 − e cos E lies in [1 − e, 1 + e], so Newton's step is always finite. The **bracket** [M − e, M + e] always contains the root because |E − M| = e|sin E| ≤ e, so bisection is always available as a fallback. From E:

- **true anomaly** ν (the angle from perigee seen from the focus): tan(ν/2) = √((1 + e)/(1 − e)) tan(E/2), computed with `atan2` so that every quadrant comes out right; ν = E for e = 0 and ν = E = π at apogee;
- **radius** r = a(1 − e cos E), from a(1 − e) at perigee to a(1 + e) at apogee;
- **position and velocity** in the perifocal frame, r = a(cos E − e, √(1 − e²) sin E), v = (√(μa)/r)(−sin E, √(1 − e²) cos E).

Three invariants check the whole chain: **vis-viva** v² = μ(2/r − 1/a), the **specific angular momentum** |r × v| = √(μa(1 − e²)), and the **specific energy** v²/2 − μ/r = −μ/2a, all constant along the orbit. For the GTO example `main()` prints: at t = 1800 s, M = 17.09°, E = 48.19°, ν = 96.88°, r = 12 551 km, v = 6869 m/s; at t = 3600 s, M = 34.18°, E = 74.36°, ν = 124.81°, r = 19 606 km, v = 4932 m/s; at apogee (t = T/2 = 18 961 s) r = 42 164 km, v = 1603 m/s; and the energy is −8.1693 MJ/kg at every time, equal to −μ/2a. The satellite covers 125° of true anomaly in the first hour and only 55° in the next 4.3 h: Kepler's second law, visible in the numbers.

### 2.6 A second root: the supersonic pitot tube

A pitot probe in supersonic flow sits behind a bow shock, and the measured ratio of stagnation pressure to free-stream static pressure is the **Rayleigh pitot formula**, a product of powers of M² with no inverse in closed form. Flight-test data reduction solves it for M at every sample. It is monotone (the ratio rises from 1.893 at M = 1 to 12.06 at M = 3 for γ = 1.4) and smooth, so a coarse bisection to 10⁻³ followed by Newton with a numerical derivative recovers M to 10⁻¹² in 17 steps; the round trip M → ratio → M is exact to 10⁻⁹ in the tests. Below the M = 1 ratio there is no shock and the formula does not apply: the solver raises rather than return a Mach number for a flow regime it does not model.

### 2.7 Stopping, and where the model stops

A root finder stops on |x_{k+1} − x_k| ≤ tol. For Newton that is safe: when the step is 10⁻¹⁴ the error is about 10⁻²⁸ (squared), so the returned iterate is good to the step size. For a linearly convergent method with rate ρ the error after stopping is about ρ/(1 − ρ) times the last step: at ρ = 0.9, nine times the tolerance. Checking |f(x)| ≤ tol instead measures something else (a residual, in the units of f) and fails badly when f′ is tiny or huge. The tests therefore check *both* the error against the oracle and the residual. Finally, the floor: no iteration can place E closer than about ulp(E) ≈ 2 × 10⁻¹⁶ rad, and the convergence-order measurement must discard errors at that floor (`P03.01`), which is what `observed_convergence_order`'s `floor` argument does. The orbit model itself stops being valid where two-body motion does: J₂ precesses the perigee of an equatorial GTO by about 0.8° per day and drag at 250 km lowers the apogee by kilometres per revolution; Kepler's equation is still solved exactly, for an orbit whose elements must be updated.

## 3. Build the computational model

Open `starter/main.py`. Every solver returns a `RootResult(root, iterates, n_evals)`; keeping the iterates is what makes the convergence measurement possible.

### Step 1: the four iterations

```python
def newton(f, dfdx, x0, tol=1e-12, max_iter=50):
    x, xs, n_evals = float(x0), [float(x0)], 0
    for _ in range(max_iter):
        fx, dfx = f(x), dfdx(x); n_evals += 2
        if dfx == 0.0:
            raise ZeroDivisionError("Newton hit a zero derivative")
        x_new = x - fx / dfx
        if not math.isfinite(x_new) or abs(x_new) > 1e150:
            raise RuntimeError("Newton diverged")
        xs.append(x_new)
        if abs(x_new - x) <= tol:
            return RootResult(x_new, np.array(xs), n_evals)
        x = x_new
    raise RuntimeError("Newton did not converge")
```

`bisection` checks the sign change (raise `ValueError`), returns an endpoint that is an exact root at once, records each midpoint and stops when the half-width is ≤ tol; `bisection_steps_needed` is (2). `secant` keeps (x_{k−1}, f_{k−1}) and (x_k, f_k), raises `ZeroDivisionError` on a flat slope and counts one evaluation per step after the first two. `fixed_point` iterates g and raises `RuntimeError` when an iterate is not finite or `max_iter` is exceeded.

### Step 2: measuring convergence

`errors_from_iterates(iterates, root)` returns |x_k − root| with exact zeros removed. `observed_convergence_order(errors, floor)` returns, for every triple of consecutive, strictly decreasing errors above the floor, q_k = log(e_{k+1}/e_k)/log(e_k/e_{k−1}). `fitted_convergence_order` fits log e_{k+1} against log e_k by least squares over all pairs above the floor. `observed_linear_rate` averages the last five ratios e_{k+1}/e_k above the floor.

### Step 3: Kepler and the orbit

`solve_kepler(M, e, method, tol)` validates e, reduces M, picks the start (E₀ = M + e sin M for Newton and fixed point, (M, M + e) for secant, [M − e, M + e] for bisection), calls the solver on `kepler_residual` and `kepler_derivative`, and shifts the result and its iterates back by 2π × revolutions. `true_anomaly` uses `2 atan2(√(1 + e) sin(E/2), √(1 − e) cos(E/2))`. `mean_motion`, `orbit_state` and `specific_energy` are the formulas of 2.5. `rayleigh_pitot_ratio` is the formula of 2.6 (raise for M < 1); `mach_from_pitot_ratio` runs `bisection` on [1, 20] to 10⁻³ then `newton` with a central-difference derivative.

Running `main()` prints the convergence table for e = 0.73, M = 60°:

```
Kepler's equation, e = 0.73, M = 60 deg: scipy brentq E = 1.763662542280359 rad
method                   E (rad)      |f(E)|  steps  f evals                         order estimates
bisection      1.763662542280358     8.9e-16     47       50  n/a: bound halves each step (rate 0.5)
fixed_point    1.763662542280359     2.2e-16     17       17                        1.00, 1.00, 1.00
secant         1.763662542280359     2.2e-16      7        7                        4.42, 1.23, 1.82
newton         1.763662542280359     2.2e-16      5       10                              2.02, 2.00
fixed-point linear rate observed 0.1399, predicted |g'(E*)| = e |cos E*| = 0.1399
bisection: 48 midpoints for tol 1e-14 on a bracket of width 1.46; predicted 48
```

then the GTO states of 2.5, the Newton step counts against eccentricity (3, 4, 7, 9, 13 steps for e = 0.1, 0.5, 0.9, 0.99, 0.999 at M = 5°), the pitot inversion, the two Newton failure demonstrations, and `root_finding.png` (error against iteration for the four methods; E(M) curves for four eccentricities).

### Verify

Run `python -m pytest curriculum/phase-03/lesson-02-root-finding`. The tests:

1. **Bisection** on cos over [0, 2]: exactly ⌈log₂(2/10⁻¹⁰)⌉ = 35 midpoints, each within 2/2ᵏ of π/2, n_evals = midpoints + 2; `ValueError` without a sign change; an exact endpoint root returned at once.
2. **Newton** on x² − 2 from 1: ≤ 6 steps, last per-step order within 0.1 of 2, equal to `scipy.optimize.newton`; raises on f′ = 0, on atan from 1.5, on the x³ − 2x + 2 cycle; linear at rate exactly ½ on (x − 1)² with more than 25 steps.
3. **Secant** on x³ − 2 equals `brentq` to 10⁻¹⁴ with one evaluation per step; raises on x² from (−1, 1).
4. **Fixed point** on cos converges to the Dottie number at observed rate |sin x*| = 0.674 ± 0.01 and diverges for g = 2x.
5. **Order estimators** recover 2, 1 and rate ½ from synthetic sequences, ignore the floor, drop zeros.
6. **Convergence order on Kepler** (e = 0.9, M = 10°): fitted orders 2.0 ± 0.2, 1.618 ± 0.2, 1.00 ± 0.05; last two Newton per-step estimates within 0.1 of 2; fixed-point rate = e|cos E*| ± 0.01; steps ordered Newton < secant < bisection < fixed point.
7. **Kepler against the oracle**: five eccentricities (0 to 0.99) × nine mean anomalies × four methods agree with `brentq` to 10⁻¹² with residuals < 10⁻¹³; e = 0 gives E = M; revolutions ±2πk are carried; e = 1, e < 0 and an unknown method raise.
8. **True anomaly**: ν = E for e = 0 modulo 2π, ν(π) = π, ν(0) = 0, ν > E near perigee, and the inverse relation cos E = (e + cos ν)/(1 + e cos ν).
9. **Orbit state** for the GTO at twelve random times over two periods: vis-viva, |r × v| and energy to 10⁻¹², r within [a(1 − e), a(1 + e)], perigee position and speed, return to perigee after T to 1 mm; n(6778 km) = 1.1314 × 10⁻³ rad/s.
10. **Rayleigh pitot**: 1.8929 at M = 1, monotone to M = 8, round trip to 10⁻⁹ at four Mach numbers, agreement with `brentq`, raises for M < 1 and for a subsonic ratio.

### Compare

`scipy.optimize.brentq` is the production scalar root finder: Brent's method combines bisection's guarantee with the secant method's speed by switching between bisection, secant and inverse quadratic interpolation steps, and it converges to full precision in about as many evaluations as the secant method while never leaving its bracket. `scipy.optimize.newton` is Newton (or the secant method when no derivative is given) with a cap and no safeguard; `root_scalar` wraps both with a common interface. On the Kepler grid `brentq` and our Newton agree to 10⁻¹² and `brentq` takes 3 to 12 evaluations against Newton's 2 to 20 (two per step). Call the library: it is tested on thousands of functions and handles the corner cases (a root at an endpoint, f returning nan) that a from-scratch routine meets only in production. What writing the four methods bought is the ability to read a convergence history and know which method produced it, to predict from f′ and f″ how many steps a solver will need, and to recognise the signatures of failure (a step count that climbs toward the cap, iterates that oscillate, an error that halves instead of squares) before they become a stale orbit in a flight computer.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Four methods, one equation.** For e = 0.73 and M = 60°, run all four methods to 10⁻¹⁴ and tabulate error against iteration (the left panel of `root_finding.png`) and against function evaluations. Predict from (2)–(5): bisection needs ⌈log₂(1.46/10⁻¹⁴)⌉ = 48 evaluations, fixed point log(10⁻¹⁴/0.084)/log(0.14) ≈ 15, Newton about 5 steps = 10 evaluations, secant about 7.
2. **The e-sweep.** For M = 5° and e = 0.1, 0.3, …, 0.9, 0.95, 0.99, 0.999, record Newton's steps from E₀ = M + e sin M and from E₀ = π, the fixed-point rate e|cos E*| and its steps, and the bisection steps. Predict where the fixed point becomes slower than bisection (rate > ½) and where Newton from E₀ = M + e sin M starts to overshoot (its first iterate exceeds π).
3. **The M-sweep.** For e = 0.9 and M from 1° to 180° in 10 values, record the same quantities. Predict that Newton is hardest near M → 0 (where f″/2f′ = e sin E/(2(1 − e cos E)) is largest) and easiest near M = 180°.
4. **Order from history.** For every (e, M) in 2 and 3, compute the fitted order of Newton and the secant method. Predict 2 and 1.618, and record the number of usable iterates each fit used; note when fewer than four errors lie above the floor and the fit cannot be made.

Raw data format: CSV tables `method, k, error_rad, evals` for 1; `e, M_deg, newton_steps_std, newton_steps_pi, fp_rate, fp_steps, bisection_steps` for 2 and 3; `e, M_deg, method, fitted_order, n_usable` for 4. Plus `root_finding.png`.

What you should find: the evaluation counts match the predictions within ±2; the fixed point crosses rate ½ near e = 0.51 at M = 5° (E* ≈ 10°, cos E* ≈ 0.98) and needs more steps than bisection beyond e ≈ 0.55 (48 against 47 at e = 0.55, 109 against 48 at e = 0.9); Newton from E₀ = M + e sin M overshoots past π for e ≥ 0.99 at M = 5° (first iterate 205°; 316° at e = 0.999) but from E₀ = π never does, at the cost of two or three extra steps at low e (5 against 3 at e = 0.1) and none at e ≥ 0.95 (8 steps for every e from 0.9 to 0.999); the fitted Newton order sits between 1.9 and 2.05 for every e ≤ 0.99 with five or more usable errors, and drops to 1.78 at e = 0.999 where the three overshooting steps pollute the fit; the secant order lies between 1.45 and 1.75; at M = 180° Newton converges in one step because E = π is an exact root and the start M + e sin M already sits on it.

## 5. Validate and improve

| quantity | theory (§2) | observed (code) | difference | cause |
|---|---|---|---|---|
| bisection midpoints, width 1.46, tol 10⁻¹⁴ | 48 | 48 | 0 | (2) is exact |
| fixed-point rate, e = 0.73, M = 60° | e|cos E*| = 0.1399 | 0.1399 | < 10⁻⁴ | (3) |
| fixed-point steps to 10⁻¹⁴ | ≈ 15 | 17 | 2 | pre-asymptotic first steps |
| Newton per-step order | 2 | 2.02, 2.00 | 0.02 | finite e_k in (4) |
| Newton steps, e = 0.99, M = 5° | — | 9 (3 wasted) | — | large f″/2f′ near perigee; overshoot to 205° |
| secant fitted order | 1.618 | 1.53 (e = 0.73), 1.61 (e = 0.9) | ≤ 0.1 | 5 to 9 usable pairs; (5) has two memory terms |
| Newton double root rate | ½ | 0.500 | < 10⁻⁹ | exact for a quadratic |
| E against `brentq`, 180 cases | 0 | ≤ 10⁻¹² | — | both at the rounding floor; `brentq`'s xtol |
| orbit invariants | exact | 10⁻¹² relative | — | rounding in √ and ÷ |
| return to perigee after T | 0 | < 1 mm | — | ulp(4 × 10⁷ m) × a few |
| pitot round trip in M | 0 | < 10⁻⁹ | — | central-difference derivative in the Newton polish |

Every row matches to within the honest tolerance, and the two that do not match *exactly* are the instructive ones: the fixed point's two extra steps and the secant method's order below 1.618 both come from the first few iterates, which no asymptotic formula describes. A convergence-order measurement is a statement about the *tail* of a sequence, and a short sequence (Newton reaches the floor in four steps) gives only one or two usable estimates. That is a limit of the measurement, not of the method, and the report should say how many iterates each estimate used.

What to improve: the starting guess. E₀ = M + e sin M is good for e ≤ 0.9 and poor near perigee for e → 1, where E₀ = π (Problem 6) or a cubic-equation estimate removes the wasted steps. What not to improve: the tolerance. Asking for 10⁻¹⁶ rad is asking for less than an ulp of E and makes every method fail its stopping test; 10⁻¹⁴ rad is 0.4 µm along a GTO at apogee, far below the orbit model's own error.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (2)–(5) and the secant order q = (1 + √5)/2; `root_finding.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a Kepler solver for a flight computer that must run at 10 Hz on orbits up to e = 0.95: method, starting guess, stopping rule, iteration cap, fallback, and the evidence for each choice.

### Adversarial self-checks

- **What if e → 1?** f′ = 1 − e cos E → 0 near perigee, Newton's constant f″/2f′ grows without bound, and the standard start overshoots; the fixed-point rate → 1 and it stalls; bisection is unaffected (48 steps whatever e). Which method would you ship for a comet-like e = 0.999 orbit, and what does that cost per solve?
- **What if the tolerance doubles, or halves?** Bisection gains or loses one step; Newton gains or loses nothing (its last step overshoots the tolerance by many orders); fixed point gains or loses log 2/log(rate) steps. A tolerance tighter than ulp(E) cannot be met by any method: check that your code fails loudly, not by looping to the cap.
- **Which assumption fails first in the real world?** Two-body motion: J₂ precesses the perigee of an equatorial GTO by about 0.8° a day, so the mean anomaly is only meaningful against regularly updated elements. The root finder is the most reliable part of the chain.
- **What measurement could falsify the model?** A ranging measurement at the predicted 22 500 km crossing time: a timing error growing linearly with time since the last element update indicates a wrong mean motion (hence a wrong a), one growing with the square indicates drag.
- **Does the result survive a change of method or start?** All four methods and `brentq` agree to 10⁻¹² from any start inside the bracket; Newton from a start outside [M − e, M + e] can still converge (the function is monotone) but may take more steps. Report E with the method and iteration count that produced it.
- **Which decision changes if the cost budget halves?** At e = 0.73 Newton uses 10 evaluations and the secant method 7 for the same tolerance; if only 5 are affordable, the secant method from a good pair still reaches 10⁻⁸ rad (0.4 m at apogee) while Newton reaches 10⁻⁶; the decision between them turns on whether f′ is cheap, and for Kepler's equation it is (one cosine).

## Check yourself

1. State the error bound of bisection after k midpoints and the number of midpoints needed to reach a tolerance tol from a bracket of width w.
2. Why does fixed-point iteration x = M + e sin x converge for every start when e < 1, and at what rate?
3. Derive the Newton step from the tangent line and state the order of its convergence at a simple root and at a double root.
4. Why is the secant method's order (1 + √5)/2, and when is it cheaper than Newton per correct digit?
5. Name three ways Newton's method can fail and give a function that exhibits each.
6. What bracket always contains the solution of Kepler's equation, and why can f′ = 1 − e cos E never vanish for e < 1?
7. A satellite has E = 90° on an orbit with e = 0.5. Give its true anomaly and its radius in units of a.
8. Which three invariants of two-body motion check a position–velocity pair produced from E, and to what precision should they hold?
