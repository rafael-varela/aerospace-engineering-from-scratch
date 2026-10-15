# P02.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: polar form, a power and three roots)

Let z = (3 + s/10) + 4i.

(a) Find |z| and arg z (in radians and degrees) by hand, then with `to_polar`.
(b) Compute z⁵ by de Moivre (modulus and argument), reduce the argument to (−π, π], and check against Python's `z ** 5`.
(c) Find the three cube roots of z with `nth_roots`; give their common modulus and their three arguments; verify that their product equals z and their sum is zero.
(d) In one sentence: why is the sum of the three roots zero?

Answer form: two numbers in (a), a complex number with its check in (b), three complex numbers and two checks in (c), one sentence in (d).

## 2. (Derivation: Euler, the addition formula, and the complex gain)

(a) Starting from the power series of e^z, cos θ and sin θ, prove e^{iθ} = cos θ + i sin θ.
(b) Using e^{i(a+b)} = e^{ia} e^{ib}, derive cos(a + b) and sin(a + b).
(c) For m ẍ + c ẋ + k x = F₀ cos ωt, substitute x = Re(X e^{iωt}) and derive X = F₀/(k − mω² + icω). State exactly where linearity and the realness of m, c, k are used.
(d) Show that |H(r)|² = 1/((1 − r²)² + 4ζ²r²) is maximal at r = √(1 − 2ζ²) for ζ < 1/√2, and give the peak value. Evaluate both for ζ = 0.08.

Answer form: four derivations; two numbers in (d).

## 3. (Unfamiliar scenario: two propellers that are *not* synchronised)

The synchrophaser is switched off and the second propeller runs slightly fast: f₁ = 40 Hz, f₂ = 40 + s/60 Hz, amplitudes 1.2 and 0.9 m/s² as in the lesson.

(a) Write x(t) = 1.2 cos(2πf₁t) + 0.9 cos(2πf₂t) as a 40 Hz signal with a slowly varying phasor: x = Re(Z(t) e^{iω₁t}). Give Z(t).
(b) Show that the envelope of x is √(1.2² + 0.9² + 2 × 1.2 × 0.9 cos(2πΔf t)) and give the beat period.
(c) Generate x(t) for 0 ≤ t ≤ 2 s at 100 kHz (synthetic) and run `fit_sinusoid` at 40 Hz over 0.1 s windows starting at t = 0, 0.5, 1.0, 1.5 s. Compare the fitted amplitudes with the envelope at the window centres.
(d) A passenger describes the cabin as "throbbing". Explain in two sentences, with your beat period, and say what the synchrophaser does to the throb.

Answer form: a formula, a period, four fitted amplitudes with four predictions, two sentences.

## 4. (Numerical, uses your code: redesign the tray mount)

The synchrophaser holds the propellers at δ = 180° − s degrees.

(a) Find the combined floor acceleration at the tray and its phase, and the corresponding floor displacement amplitude at 40 Hz.
(b) With the existing 30 Hz, ζ = 0.08 mount, find the transmissibility and the tray acceleration.
(c) Use `mount_frequency_for_transmissibility` to find the mount frequency for T = 0.2 at 40 Hz. Give the new k and c for the 2.5 kg tray and the static sag of the tray under gravity.
(d) Check the new mount at 80 Hz (blade passing of a two-blade propeller) and at 20 Hz (1200 rpm idle). Then verify the 40 Hz result by integrating the base-excited ODE m ẍ + c ẋ + k x = c ẏ + k y with `solve_ivp` for 3 s and fitting the last second.

Answer form: three numbers in (a), two in (b), four in (c), three transmissibilities and one fitted ratio in (d).

## 5. (Judgment: can one synchrophaser angle quiet the whole cabin?)

At the tray the propeller phasors are 1.2 at 0° and 0.9 at δ. Two rows aft, the structural paths are different: the first propeller arrives at 1.0 m/s² and 0° and the second at 1.0 m/s² and δ + 60°.

(a) For δ = 180°, compute the combined amplitude at both locations.
(b) Find the δ that minimises each location separately, and the δ that minimises the *sum* of the two amplitudes (a coarse sweep with `combined_amplitude` is enough).
(c) The phasor model assumes the structure is linear and that the 40 Hz line is all that matters. Name two things a real cabin measurement would contain that the model does not, and say for each whether it would change the recommended δ.
(d) Would you recommend a fixed δ or a δ schedule? Justify with your numbers.

Answer form: two amplitudes, three angles, two named effects, one recommendation.

## 6. (Hard: why a balanced propeller has no once-per-revolution force)

Each blade of an N-bladed propeller produces, at the hub, a force whose once-per-revolution component has amplitude F and phase 2πk/N for blade k (the blades are identical and equally spaced).

(a) Show that the sum Σ_{k=0}^{N−1} e^{2πik/N} is zero for N ≥ 2, so the 1/rev hub force vanishes.
(b) The harmonic h of blade k has phase 2πhk/N. Using the geometric-series formula |Σ_{k=0}^{N−1} e^{ikδ}| = |sin(Nδ/2)/sin(δ/2)|, show that the hub force at harmonic h is zero unless h is a multiple of N, where it is NF. Confirm with `add_sinusoids` for N = 2, 3, 4 and h = 1, …, 6.
(c) One blade of a three-blade propeller is 5 % heavier (its force is 1.05F at phase 0). Compute the residual 1/rev hub force with `add_sinusoids`. Generalise: the residual equals the imbalance of one blade.
(d) Derive the geometric-series formula in (b) from the sum of a geometric progression and Euler's formula.

Answer form: two proofs, a table of six harmonics for three blade counts, one residual, one derivation.
