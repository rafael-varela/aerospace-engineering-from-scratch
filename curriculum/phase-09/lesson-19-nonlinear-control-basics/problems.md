# P09.19  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a Lyapunov equation by hand)

A damped oscillator ẋ = Ax with A = [[0, 1], [−k, −c]], k = 4 + s/5 s⁻², c = 1 s⁻¹.

(a) Solve AᵀP + PA + I = 0 by hand: write the three scalar equations for p₁₁, p₁₂, p₂₂ and solve them. (b) Check that P is positive definite and that V = xᵀPx decreases along trajectories at the rate xᵀx. (c) Check (a) with `lyapunov_solve` and `scipy.linalg.solve_continuous_lyapunov`. (d) What happens to P as c → 0, and why does that make sense?

Answer form: three entries with a derivation, two checks, two agreeing matrices, one limit with one sentence.

## 2. (Derivation: the four tools)

(a) Show that if A is Hurwitz, P = ∫₀^∞ e^{Aᵀt}Qe^{At} dt solves AᵀP + PA + Q = 0 and is positive definite for Q > 0; and conversely that P > 0 with Q > 0 implies A is Hurwitz. (b) Derive the describing function of a saturation of slope k and linear range δ, N(A) = (2k/π)[arcsin(δ/A) + (δ/A)√(1 − (δ/A)²)] for A > δ, from the first Fourier coefficient. (c) Derive the describing function of a relay ±M with hysteresis h, and show that −1/N(A) lies on the horizontal line Im = −πh/(4M). (d) For m l²θ̈ = mgl sin θ − bθ̇ + u, derive the feedback-linearising law that gives θ̈ = v, and state precisely what it needs to know.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a thruster limit cycle)

A sounding rocket's roll is controlled by a pair of cold-gas thrusters through an on-off valve with hysteresis. From valve command to roll angle the plant is G(s) = K/(s(τs + 1)) with K = 2 + s/10 rad/s² per unit command and valve lag τ = 0.05 s; the valve is a relay ±1 with hysteresis h = 0.01 rad on the roll error.

(a) Show that without hysteresis harmonic balance predicts no limit cycle (the phase of G never reaches −180°), yet the loop does not settle to zero either. (b) With h = 0.01 rad, predict the limit-cycle amplitude and frequency with `predict_relay_hysteresis_limit_cycle`. (c) Simulate with `simulate_relay_loop` and measure both. (d) How would you reduce the amplitude, and what does it cost in thruster firings per second?

Answer form: one argument, a predicted (A, ω), a simulated (A, ω), one design change with a cost.

## 4. (Numerical, uses your code: your loops)

(a) For the fin servo of the lesson with a saturating amplifier of limit 1 V and gain m × k_crit, m = 1.5 + s/20, predict the limit cycle and measure it in simulation. (b) For the pendulum, repeat the region-of-attraction estimate with Q = diag(1, s/100) and Q = diag(1, 1). Report c*, the angle reached at rest inside each estimate, and the fraction of 100 seeded random starts inside each estimate that converge.

Answer form: a predicted and simulated (A, ω) pair with percentage errors; a two-row table.

## 5. (Judgment: trusting nonlinear tools)

(a) The describing function predicted the relay limit cycle within 5 %. Name two plants for which you would not trust it, and why. (b) The region-of-attraction estimate reached ±23.6° at rest while simulation recovers up to 37.7°. Is the estimate wrong? What would make it less conservative? (c) The feedback-linearising controller from 60° needs 16.9 N·m, but the stand's motor gives 6 N·m. What happens, and what does that say about feedback linearisation in practice? (d) A pendulum whose mass is 10 % larger than modelled is feedback-linearised with the nominal mass. Is the closed loop still stable near upright? Give a reason, not a simulation.

Answer form: four short paragraphs.

## 6. (Hard: feedback linearisation with the wrong model and a torque limit)

(a) Simulate the feedback-linearised pendulum from 60° (the lesson's ω_n = 4 rad/s, ζ = 0.7) when the true mass is 1.1 and 1.2 times the modelled one, and report the largest deviation from the linear closed form. (b) Add the 6 N·m torque limit to both the feedback-linearising and the linear controllers (`simulate_pendulum` with a clipped controller) and find, by bisection, the largest initial angle at rest from which each recovers. (c) Explain the result of (b) using m g l sin θ and the torque limit.

Answer form: two deviations, two angles, one explanation.
