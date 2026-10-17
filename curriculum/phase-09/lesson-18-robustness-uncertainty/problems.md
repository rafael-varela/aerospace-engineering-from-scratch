# P09.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: peaks and the margins they guarantee)

(a) A loop has Ms = 1.5 + s/30. Give the disk-margin bounds: the gain range (Ms/(Ms + 1), Ms/(Ms − 1)) and the phase margin 2 arcsin(1/(2Ms)). (b) Give ‖G‖∞ by hand for G₁ = 3/(s + 2) and for G₂ = ω_n²/(s² + 2ζω_ns + ω_n²) with ω_n = 5 rad/s and ζ = 0.1 + s/100, and the frequency of each peak. (c) Check (b) with `tf_hinf_norm` by both methods. (d) A designer says "Ms = 2 means a 6 dB gain margin". Is that the guarantee? State it correctly in dB.

Answer form: two bounds and an angle; two norms with frequencies; four agreeing numbers; one corrected sentence.

## 2. (Derivation: the four results the lesson rests on)

(a) Show S + T = 1. For a reference r, an output disturbance d (added to y) and sensor noise n (added to the measurement), write the output y, the tracking error r − y and the control u in terms of S, T and KS. (b) Show that min_ω |1 + L(jω)| = 1/Ms, and derive the gain and phase guarantees of the disk margin from the geometry of the disk of radius 1/Ms about −1. (c) For G = G₀(1 + WΔ), show 1 + CG = (1 + L₀)(1 + WTΔ) and deduce, with the Nyquist criterion and |WTΔ| < 1, that ‖WT‖∞ < 1 guarantees stability for every stable Δ with ‖Δ‖∞ ≤ 1. Then show it is also necessary by constructing a destabilising all-pass Δ when ‖WT‖∞ ≥ 1. (d) For G = C(sI − A)⁻¹B with D = 0, show that jω is an eigenvalue of H(γ) = [[A, BBᵀ/γ²], [−CᵀC, −Aᵀ]] if and only if γ is a singular value of G(jω).

Answer form: four derivations, each ending in a boxed result.

## 3. (Unfamiliar scenario: the waterbed of an unstable fighter)

A relaxed-stability aircraft's pitch mode diverges at p = 1 + s/10 rad/s without feedback. The designer wants |S(jω)| ≤ 0.3 (disturbance rejection) for ω ≤ 2 rad/s, and the actuator and structure make |S| ≈ 1 above ω₂ = 30 rad/s whatever the controller does.

(a) Use the Bode integral ∫₀^∞ ln|S| dω = πp to show that ln Ms ≥ (πp − 2 ln 0.3)/(30 − 2), and evaluate the bound on Ms. (b) What does it become if the actuator bandwidth halves (ω₂ = 15 rad/s)? (c) Check the integral numerically for the lesson's `unstable_airframe` with two stabilising lead controllers of your choice (`bode_sensitivity_integral`). (d) Explain in two sentences why no controller can remove the bound, and what the airframe designer (not the control designer) could change.

Answer form: a bound with its arithmetic, a second bound, two integrals against πp, two sentences.

## 4. (Numerical, uses your code: your servo design)

Use the lead C(s) = K_c(s/20 + 1)/(s/80 + 1) with K_c = 0.4 + s/50 on the nominal servo.

(a) Give Ms, Mt, ‖KS‖∞ and the disk-margin bounds. (b) Give ‖S‖∞ by `hinf_norm_sweep` and by `hinf_norm_bisection`, to eight digits. (c) Give ‖WT‖∞ and say whether robust stability is guaranteed for the lesson's weight. (d) Measure the phase and gain margins of the *true* servo loop (exact delay) and check them against the disk bounds computed from the true loop's own Ms.

Answer form: four numbers and a range; two agreeing numbers; one number with a verdict; two margins with two bounds.

## 5. (Judgment: what the weight means)

(a) The design with K_c = 1.6 fails the small-gain test (‖WT‖∞ = 1.108) yet the true servo loop is stable with a 39° phase margin. Is the test wrong? Explain what it guarantees and what it does not. (b) The weight was fitted to one actuator. A second unit off the production line has its compliance mode at 120 rad/s instead of 150 rad/s. Is it covered? Check with `weight_covers` (modify `true_servo_relative_error` in a copy). (c) The real uncertainty in the motor gain is a real number between 0.9 and 1.1, but the small-gain test allows any complex Δ. Which way does that make the test err, and by roughly how much for the chosen design? (d) Would you fly the K_c = 1.2 design (‖WT‖∞ = 0.974)? Name what you would measure first.

Answer form: four short paragraphs with specific numbers.

## 6. (Hard: the robust gain limit against the true limit)

(a) With z = 20 rad/s and p = 80 rad/s fixed, find by bisection the largest K_c for which ‖WT‖∞ < 1 (the robustly guaranteed gain). (b) Find by bisection the largest K_c for which the true servo (Padé order 4 for the delay) is closed-loop stable. (c) For K_c 5 % above the value in (a), construct the destabilising Δ with `destabilising_perturbation` and report its size ‖Δ‖∞ and the frequency of the resulting closed-loop pole on the imaginary axis. (d) Explain the gap between (a) and (b) in terms of the weight and the phase of the true Δ.

Answer form: two gains, one size and one frequency, a short explanation.
