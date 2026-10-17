# P06.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, report the energy residual of every simulation you use, and say which numbers are illustrative.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the gas spring)

Use the lesson's strut (A_a = 7.854 × 10⁻³ m², p₀ = 1.5 MPa, V₀ = 3.1416 × 10⁻³ m³).

(a) Compute the gas force at s = 0, 0.15 and 0.30 m for n = 1.0 (isothermal; use 1.0001 in the code), 1.1 and 1.4, and the stored energy at 0.30 m for each. Check the isothermal energy with the closed form p₀V₀ ln(V₀/V) − p_atm A_a s.
(b) An airframe share of M = 1100 + 20s kg sits on the strut at rest on the ground with no lift. Find the static stroke (gas force = Mg cos β) and the gas spring's stiffness dF/ds there.
(c) What natural frequency would the airframe have bouncing on that gas spring alone? Why is a landing gear so softly sprung compared with the isolator of `P06.17`?

Answer form: nine forces and three energies with a check; a stroke and a stiffness; a frequency and a sentence.

## 2. (Derivation)

(a) For the reduced coordinates (y₁, s), derive the axle kinematics, the kinetic energy and the mass matrix M_r, and the generalised forces Q_y₁ and Q_s by virtual work. Explain why the strut force appears only in Q_s and why the joint forces appear in neither.
(b) Derive the gas energy U_gas(s) in closed form, and its isothermal limit.
(c) For the full coordinates, write the four constraints, derive the Jacobian Φ_q and the acceleration right-hand side γ₄ for Φ₄ = (r₂ − r₁)·n(θ₁).
(d) Prove that the constraint forces Φ_qᵀλ do no work on any motion that satisfies the constraints. What does that imply for the energy balance of the Newton–Euler model?

Answer form: four derivations.

## 3. (Unfamiliar scenario: a crushable lander leg)

A 300 kg Mars lander touches down at 2.5 m/s on four legs with crushable aluminium honeycomb absorbers (a constant crushing force, no rebound) and 0.15 m of crush stroke each. Mars gravity is 3.711 m/s²; there is no lift; take the four legs as equally loaded and the motion as vertical.

(a) How much energy must the legs absorb, including the work of gravity over the stroke? What constant force per leg uses exactly the full stroke?
(b) What is the deceleration, in Earth g, and the load factor (total leg force divided by the lander's Earth weight)?
(c) What is the efficiency of a crushable absorber, and how much larger would the peak force be with a linear spring of the same stroke? Why do planetary landers use crushables while aircraft use oleos?

Answer form: an energy and a force; two numbers; an efficiency, a ratio and a sentence.

## 4. (Numerical, uses your code: choose the orifice)

For a sink rate of v = 2.5 + s/20 m/s with lift = weight:

(a) Sweep the orifice area from 3.0 to 6.0 cm² in steps of 0.1 cm² and find the area that minimises the gear load factor. Report the load factor, stroke, efficiency and maximum energy residual at the optimum, and the same quantities for the lesson's 4.0 cm² orifice.
(b) Plot the load factor against orifice area and explain the shape (what dominates on each side of the minimum?).
(c) Would you choose the exact optimum? Consider the C_d uncertainty of experiment 3 in the lesson.

Answer form: two rows of four numbers; a plot and an explanation; a decision.

## 5. (Numerical, uses your code: when does it bottom?)

(a) Find, with `scipy.optimize.brentq` on the peak stroke minus 0.30 m, the sink rate at which the lesson's strut just bottoms, with lift = weight and with no lift. Report the load factor at each.
(b) Explain why the two limits differ by the amount they do, using the energy balance.
(c) The certification requirement in the lesson is 3.05 m/s with lift = weight. What margin in sink rate does the strut have? Is that the right way to express the margin?

Answer form: two sink rates with load factors; an energy argument; a margin and a sentence.

## 6. (Judgment: the drop test disagrees)

A drop test of the lesson's gear at 3.05 m/s with lift = weight records a peak ground force of 36.5 kN at 16 ms after touchdown, against the model's 32.2 kN at 34 ms.

(a) Run the model with C_d = 0.6 and with the tyre stiffness raised to 1.2 × 10⁶ and 1.6 × 10⁶ N/m. Report the peak force and its time for each.
(b) Which explanation does the *timing* support? Which does it rule out, even though it can match the peak force?
(c) Name two further measurements from the same drop test that would confirm your explanation, and say what each should show.
(d) Should the model's C_d be tuned until the peak force matches? Argue.

Answer form: three pairs of numbers; an argument; two measurements; a decision.
