# P02.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: classify by hand)

A = [[−1, 2], [−s/10, −3]] (units s⁻¹).

(a) Find τ and Δ and classify the fixed point from the trace–determinant plane. (b) Find the eigenvalues. If they are complex, give ω_n, ζ, the period and the time to half amplitude; if real, give both time constants and the eigenvectors. (c) Give the slopes of the two nullclines and state the sense of rotation (clockwise or anticlockwise) of trajectories near the origin. (d) Confirm with `classify`, `eigen_structure`, `modal_parameters` and `nullcline_slopes`.

Answer form: two numbers, a type, two eigenvalues, four modal numbers or two eigenvectors, two slopes, one sense.

## 2. (Derivation: the plane, the real form, and the phugoid)

(a) Starting from **x** = **v**e^{λt}, derive det(A − λI) = λ² − τλ + Δ for a 2×2 A and show λ₁ + λ₂ = τ, λ₁λ₂ = Δ.
(b) For a complex pair λ = σ ± iω with eigenvector **v** = **p** + i**q**, show that e^{σt}(**p** cos ωt − **q** sin ωt) and e^{σt}(**p** sin ωt + **q** cos ωt) are real solutions and that the motion is an ellipse in the (**p**, **q**) frame scaled by e^{σt}.
(c) Prove that for τ = 0 the quadratic Q = c x₁² − 2a x₁x₂ − b x₂² is constant along trajectories, and show that Q is a positive- or negative-definite form exactly when Δ > 0 (a center) and indefinite when Δ < 0 (a saddle: the level curves are hyperbolas).
(d) Derive Lanchester's phugoid matrix from the two statements "excess speed u raises lift by the fraction 2u/V₀, which curves the path" and "a climb angle θ decelerates the aircraft at gθ, and excess drag removes energy at 2g u/(V₀ L/D)", and show ω_n = √2 g/V₀, ζ = 1/(√2 L/D).

Answer form: four derivations.

## 3. (Unfamiliar scenario: propellant transfer between two tanks)

Two tanks exchange propellant through a cross-feed line: ṁ₁ = −k₁m₁ + k₂m₂, ṁ₂ = k₁m₁ − k₂m₂ with k₁ = 0.010 s⁻¹ and k₂ = 0.005 + s/2000 s⁻¹. Initially m₁ = 100 kg, m₂ = 0.

(a) Find τ and Δ and classify. Why is this fixed point not isolated, and what does the line of fixed points mean physically? (b) Find the eigenvalues and eigenvectors by hand. Which eigenvector is the line of equilibria and which is the direction of approach? (c) Write m₁(t), m₂(t) in closed form, give the final masses and the time constant, and confirm with `solve_linear` at t = 0, 100, …, 600 s. (d) Which quantity is conserved exactly along every trajectory, and how does it show up in the eigenstructure?

Answer form: two numbers, a type with an explanation, two eigenpairs, two formulas with a time constant, one conserved quantity.

## 4. (Numerical, uses your code: your aircraft's pitch damping)

Take the short-period model with Z_α/V₀ = −0.8 s⁻¹, M_α = −4.0 s⁻², and M_q = −(0.6 + s/25) s⁻¹.

(a) Classify and give ω_n, ζ, period and time to half amplitude. (b) Propagate a 5° angle-of-attack disturbance (q(0) = 0) for 8 s with `solve_linear` and report the peak pitch rate (deg/s) and its time, and the time after which |α| stays below 0.5°. (c) Repeat with `rk4_trajectory` at 800 steps and give the maximum difference. (d) Draw the phase portrait from eight starts on a circle of radius 0.3 and the two nullclines; mark where trajectories cross each nullcline vertically and horizontally.

Answer form: a type, four numbers, two times and a rate, one difference, one figure.

## 5. (Judgment: what the 2×2 models leave out)

(a) The short-period model freezes speed and the phugoid model freezes angle of attack. Using the two periods (3.2 s and 104 s for the lesson's aircraft), argue why each approximation is good here and name an aircraft type for which it would not be. (b) The saddle at M_α = +2 doubles every 1.62 s. What does a pilot experience, and why is "unstable" not the same as "unflyable" for a modern airliner? (c) A flight-test pulse gives a short-period damping ratio of 0.35 where the model predicts 0.45. Name two physical effects the model omits that could account for the difference and say which way each pushes ζ. (d) State one measurable sign that the linear model has stopped applying during a large disturbance.

Answer form: four short paragraphs with specific mechanisms.

## 6. (Hard: walking into the saddle, and walking back)

With Z_α/V₀ = −0.8 and M_q = −1.2 fixed, the short-period matrix depends on M_α alone.

(a) Derive the two transition values of M_α (spiral → node, node → saddle) and verify them with `classify` on a fine sweep. (b) For M_α = 1.0, 1.5, 2.0 and 3.0 give the unstable eigenvalue, the doubling time and the angle of the unstable eigenvector in the (α, q) plane; explain the trend of the angle. (c) A stability augmentation system adds pitch damping by feeding back q to the elevator, changing M_q to M_q + k with k < 0. For M_α = +2, show algebraically that *no* value of k gives ζ = 0.5 (the fixed point becomes a node, never a well-damped spiral), and explain why in terms of Δ. (d) The alternative is to feed back α (changing M_α). What effective M_α restores ζ = 0.45 with M_q = −1.2, and what does that say about what an augmentation system for a relaxed-stability aircraft must sense? Verify with `modal_parameters`.

Answer form: two derivations with checks, a four-row table with an explanation, an algebraic argument, one number with a conclusion.
