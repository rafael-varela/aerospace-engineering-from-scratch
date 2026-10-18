# P10.19  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed form) and which are numerical. Use Theodorsen's conventions from the lesson: b is the semi-chord, the pivot is a semi-chords aft of mid-chord, plunge h is positive down, s = Vt/b.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: how unsteady is it?)

Compute the reduced frequency k = ωb/V for each case, say whether the flow is quasi-steady (k < 0.05), unsteady (0.05 ≤ k ≤ 0.2) or highly unsteady (k > 0.2), and give |C(k)| and its phase from `theodorsen`:

(a) a sailplane wing, chord 0.9 m at 25 m/s, with a 3 Hz torsional vibration;
(b) the reference trainer's wing (chord 1.2 m) at 55 m/s with a 5 Hz torsional vibration;
(c) a section of the reference helicopter rotor (chord 0.25 m, tip speed 210 m/s, radius 5.0 m) at 0.75 R, whose pitch varies once per revolution (cyclic pitch).

For (c), also give the amplitude of the total lift relative to 2πα̂ and its phase for pitching about the quarter chord.

Answer form: three k values with a regime each, three |C| with phases; one ratio and one phase.

## 2. (Derivation: from the indicial function to C(k))

(a) Write Duhamel's integral for the circulatory lift of a section whose three-quarter-chord angle α₃/₄(s) starts from zero, using the Wagner function φ(s). Substitute α₃/₄ = e^{iks} (started at s = −∞) and show that the steady-state lift is 2π α₃/₄ · ik Φ(ik), where Φ(p) is the Laplace transform of φ (`P02.08`).
(b) For φ(s) = 1 − Σ Aᵢ e^{−βᵢs}, find Φ(p) and show that C(k) = 1 − Σ Aᵢ ik/(ik + βᵢ). Find its limits as k → 0 and k → ∞ and say what each means physically.
(c) Integrate by parts to show that C(k) = φ(0) + ∫₀^∞ φ′(τ) e^{−ikτ} dτ. Why is this form better for a numerical check than ik∫₀^∞ φ e^{−ikτ} dτ?
(d) For a one-term indicial function φ = 1 − ½e^{−βs}, find the reduced time at which the lift reaches 90 % of its steady value, and choose β so that it matches Jones's 12.7 semi-chords. What does a one-term fit get wrong?

Answer form: four short derivations; two limits; a value of β with a sentence.

## 3. (Unfamiliar scenario: a rotor blade in cyclic pitch)

A helicopter in forward flight varies blade pitch once per revolution. Treat a section at 0.75 R of the reference rotor (Problem 1(c)) as a two-dimensional section pitching about its quarter chord with an amplitude of 3°, in a steady stream of 157.5 m/s (forward speed ignored).

(a) Compute the lift-coefficient amplitude by quasi-steady theory and by Theodorsen. By how many degrees of azimuth does the lift lag or lead the pitch?
(b) The blade's flapping response to that lift is itself a lagging second-order system. Explain why a few degrees of aerodynamic phase matter to the helicopter's control designer even though the amplitude change is only a few per cent.
(c) Name two effects that the two-dimensional, attached-flow Theodorsen model leaves out for a real rotor, and say which way each changes the answer.

Answer form: two amplitudes, a phase in degrees of azimuth; two short paragraphs.

## 4. (Numerical, uses your code: a pitching section of your own)

Let k = 0.1 + s/100 and the pivot a = −0.5 + s/100 (so the pivot moves aft of the quarter chord as s grows).

(a) Compute the complex lift per unit α̂ from Theodorsen and from quasi-steady theory. Give both amplitudes and phases.
(b) Simulate the motion from rest with `simulate_harmonic_pitch` (α̂ = 2°) and show that after the transient its amplitude and phase agree with Theodorsen to within the error of Jones's approximation. State that error at your k.
(c) Plot the lift loop for α = 4° ± 2°. Compute ∮ cl dα from the loop and from π Im(ᾱ* ĉl), and say what the sign of the loop's direction tells you.
(d) Compute the aerodynamic work per cycle on the pitching motion (the moment about your pivot). Is the motion damped by the air?

Answer form: four complex numbers as amplitude and phase; an agreement statement with an error; two equal numbers and a sentence; a number and a yes/no.

## 5. (Numerical, uses your code: gust alleviation)

The trainer (chord 1.2 m, 55 m/s) flies through a vertical "1 − cos" gust of peak 5 m/s with gradient length H = (5 + s/5) chords (gust length 2H).

(a) Compute the peak section lift increment by quasi-steady theory and with Küssner's function (`gust_response`). Report the ratio.
(b) Repeat for H = 12.5 chords. Why does the unsteady reduction shrink as the gust gets longer?
(c) The trainer's wing loading is 750 × 9.81/12.0 N/m². Using the finite-wing lift slope 4.63 per radian (from `P10.20`'s lattice, extrapolated to convergence) in place of 2π, estimate the peak load-factor increment Δn for your gust, with and without the section's unsteady reduction. What else (that this lesson leaves out) reduces Δn further in a real aircraft?

Answer form: two lift increments and a ratio; a ratio and a sentence; two Δn values and a sentence.

## 6. (Judgment: where does single-degree-of-freedom pitch flutter come from?)

(a) With `pitch_flutter_boundary`, find the reduced frequency below which a section pitching about a = −0.6, −0.7, −0.8, −1.0 is fed energy by the air. Express each pivot as a fraction of the chord from the leading edge.
(b) Explain, using the split of the moment into an apparent-mass part and a circulatory part acting at the quarter chord, why pivots at and behind the quarter chord are always damped.
(c) Real flutter of an attached-flow wing usually couples bending (plunge) and torsion (pitch). Using Problem 4's result that pure plunge is always damped, argue why neither single-degree-of-freedom result predicts classical flutter, and what extra information (structural, not aerodynamic) a flutter calculation needs (`P13.18`).
(d) Which assumption of Theodorsen's theory fails first for a helicopter blade or a wind-turbine blade near stall, and what happens to the damping then?

Answer form: four k values and four chord fractions; three short arguments.
