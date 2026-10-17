# P03.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are bounds and which are observed.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: two Newton steps by hand)

A Molniya-type orbit has e = 0.74. The mean anomaly is M = (10 s) mod 360 degrees.

(a) Convert M to radians, form E₀ = M + e sin M, and carry out two Newton steps by hand (calculator allowed), writing f, f′ and the step at each stage.
(b) Give the residual f(E₂) and estimate the error of E₂ from the quadratic-convergence rule, using the size of the second step.
(c) Check with `solve_kepler` and report how many steps the code took to 10⁻¹⁴ and the true anomaly in degrees.
(d) How many bisection midpoints would the bracket [M − e, M + e] need for the same tolerance?

Answer form: two hand steps with intermediate values, a residual, an error estimate, three numbers from the code.

## 2. (Derivation: the four error recurrences)

(a) From f(x*) = 0 and a second-order Taylor expansion about x_k, derive e_{k+1} = |f″(ξ)/(2f′(x_k))| e_k² for Newton's method, and state what is assumed about f′(x*).
(b) Show that at a double root, f(x) = (x − x*)² h(x) with h(x*) ≠ 0, Newton's recurrence becomes e_{k+1} = e_k/2 + O(e_k²), and that for f = (x − x*)² exactly it is e_{k+1} = e_k/2 with no remainder.
(c) Derive e_{k+1} ≈ |g′(x*)| e_k for fixed-point iteration from the mean value theorem, and give the condition for convergence.
(d) For the secant method, show that e_{k+1} ≈ |f″/2f′| e_k e_{k−1} (use the Newton divided-difference form of the interpolation error through two points), then assume e_{k+1} = C e_k^q and derive q² = q + 1.

Answer form: four derivations; the final value of q.

## 3. (Unfamiliar scenario: a hyperbolic flyby)

A spacecraft on a hyperbolic flyby has eccentricity e = 1.2 + s/100 > 1. The hyperbolic Kepler equation is M_h = e sinh H − H, where H is the hyperbolic anomaly and M_h = 2.0 at the time of interest.

(a) Show that f(H) = e sinh H − H − M_h is monotone increasing for e > 1 and has exactly one root; give a bracket that must contain it.
(b) Solve for H with your `newton` from H₀ = asinh(M_h/e) and from H₀ = M_h; report H, the step counts, and check with `scipy.optimize.brentq`.
(c) The true anomaly follows from tan(ν/2) = √((e + 1)/(e − 1)) tanh(H/2). Give ν in degrees and compare with the asymptotic true anomaly ν_∞ = arccos(−1/e) that the spacecraft approaches as t → ∞.

Answer form: a short argument, two numbers with step counts, two angles.

## 4. (Numerical, uses your code: measured orders)

Take e = 0.5 + s/100 and M = 20°.

(a) Solve with all four methods to 10⁻¹⁴ and tabulate steps and function evaluations.
(b) For each method compute the errors against `scipy.optimize.brentq` (xtol 10⁻¹⁵), then the fitted convergence order and, for fixed point, the observed linear rate. Compare with theory (2, 1.618, 1 with rate e|cos E*|, and bisection's bound).
(c) Which method is cheapest per correct digit here? Would the answer change if f′ cost ten times more than f?

Answer form: a 4-row table (steps, evaluations, fitted order, theory), one paragraph.

## 5. (Numerical, uses your code: when does the GTO satellite cross a given radius?)

The GTO of the lesson has perigee 250 km, apogee 35 786 km. Find the first time after perigee at which the radius equals R = (15 000 + 500 s) km.

(a) Write g(t) = |r(t)| − R using `orbit_state` and solve g(t) = 0 by bisection on [0, T/2] to 10⁻⁶ s and by the secant method from (1000 s, 2000 s). Report both times, their difference, and the number of Kepler solves each needed.
(b) Check analytically: from r = a(1 − e cos E) find E at the crossing, then M = E − e sin E and t = M/n. How close is the root finder to this?
(c) When is the second crossing, and why does bisection on [0, T] fail to find a unique one?

Answer form: three times in seconds and minutes, a solve count, one explanation.

## 6. (Judgment: specifying the flight-computer solver)

A guidance computer must solve Kepler's equation at 10 Hz for orbits with e up to 0.95, within a budget of 12 evaluations of sin/cos per solve, and must never return a wrong E.

(a) Which method and starting guess do you choose, and why not each of the other three?
(b) Specify the stopping rule, the iteration cap and what happens when the cap is hit.
(c) The requirement is extended to e up to 0.999 for a comet-rendezvous mission. What changes, and what does it cost?
(d) A reviewer asks for |f(E)| < 10⁻¹² as the stopping criterion instead of the step size. Is that better, worse or equivalent here? Justify with f′'s range.

Answer form: four short paragraphs, each with at least one number.

## 7. (Proof: the Kepler map is a contraction)

(a) Prove that g(E) = M + e sin E satisfies |g(E₁) − g(E₂)| ≤ e |E₁ − E₂| for all E₁, E₂ when 0 ≤ e < 1.
(b) Deduce that Kepler's equation has exactly one solution for every M, and that fixed-point iteration converges to it from any starting value, with |E_k − E*| ≤ eᵏ |E₀ − E*|.
(c) Use (b) to bound the number of iterations needed to reach 10⁻¹² from E₀ = M + e sin M for e = 0.74, and compare with the count your code produces for the orbit of Problem 1. Why is the bound pessimistic?

Answer form: a proof, a deduction, one number compared with one observation.
