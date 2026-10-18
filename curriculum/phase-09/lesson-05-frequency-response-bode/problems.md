# P09.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Every sine test you run is a SIMULATION on a model: label it so in every table and plot.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a servo loop by hand)

L(s) = K/(s(τs + 1)) with K = s s⁻¹ and τ = 0.1 s.

(a) Write L in Bode form and give K_B, the type and the corner. (b) Sketch the straight-line magnitude and phase, and read off their values at ω = 1, 10 and 100 rad/s. (c) Compute the exact values with `bode` and explain each difference. (d) Find the gain crossover ω_c in closed form and with `gain_crossover`, and the phase there. Where does the straight line cross 0 dB?

Answer form: a Bode form, a 3 × 4 table, a frequency and a phase, one straight-line crossover.

## 2. (Derivation: the frequency response)

(a) Derive the steady-state response (1) from the convolution integral, and show that the transient decays at the rate of the slowest pole. (b) Show that dB and phase add in a cascade. (c) For a real pole, derive the −3.01 dB corner error, the −0.043 dB error a decade away and the 5.71° maximum error of the decade phase rule. (d) Derive ω_r and M_r of (3) and the magnitude error −20 log₁₀(2ζ) of the straight line at ω_n.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a notch filter for a wing-bending mode)

A rate gyro on a small aircraft picks up the first wing-bending mode at f_m = 8 + s/5 Hz with damping ratio 0.02, modelled as ω_m²/(s² + 2·0.02 ω_m s + ω_m²).

(a) Give the resonant peak of the mode in dB. (b) A notch filter (s² + ω_m²)/(s² + 2ζ_n ω_m s + ω_m²) with ζ_n = 0.5 is placed after the gyro. Compute its magnitude and phase at 1 Hz, 2 Hz, f_m ± 0.5 Hz and 2f_m. (c) Find the largest magnitude of mode × notch over 1–40 Hz and its frequency. (d) What does the notch cost the pitch loop, whose crossover is near 1 Hz, and what happens if the mode frequency is 10 % lower than modelled (fuel burned, a different store)?

Answer form: a number, a five-row table, a peak with its frequency, a paragraph.

## 4. (Numerical, uses your code: how long must a sine test wait?)

Use the closed pitch-attitude loop of `P09.03`, θ/θ_c = 720(s + 1.5)(s + 100)/(s⁵ + 123.2s⁴ + 2400s³ + 16 240s² + 115 880s + 108 000), at ω = 1 + s/10 rad/s.

(a) Compute the exact |G| (dB) and ∠G. (b) Run `sine_test` with settle = 0, 1, 3 and 10 s and three fitted periods, and tabulate the gain and phase errors. (c) Explain the trend from the pole locations. (d) How long would you wait at each frequency of a 0.1–10 rad/s sweep to keep the phase error below 0.01°, and how long does the whole sweep take with 20 frequencies?

Answer form: two numbers, a four-row table, one sentence, a time per point and a total.

## 5. (Judgment: reading and trusting a Bode plot)

(a) The straight-line Bode of the pitch loop is off by up to 4.2 dB and 22.6°. For which decisions is it good enough, and for which not? (b) Two plants have identical magnitude plots, one with a zero at −2 and one at +2. Which phase would you expect from each, and why does the difference matter for a controller designed from magnitude only? (c) A ground test of the elevator servo uses a ±0.2 deg sine and a technician reports "the gain is 3 dB below the model at every frequency". Name two physical causes, and the test you would run to decide between them. (d) Why must the open-loop frequency response of a statically unstable airframe be obtained from a closed-loop test?

Answer form: four short paragraphs.

## 6. (Hard: a time delay)

A flight computer adds a computation and transport delay T to the servo loop L(s) = 14.31/(s(0.05s + 1)), so L_T(s) = L(s)e^{−sT}.

(a) Show that e^{−jωT} has magnitude 1 and phase −ωT, and that the delay leaves the gain crossover unchanged. (b) Compute ω_c and the phase there without delay, and the phase at ω_c for T = 10, 20 and 50 ms. What is 180° + phase in each case? (c) The first-order Padé approximation e^{−sT} ≈ (1 − sT/2)/(1 + sT/2) is an all-pass with a right-half-plane zero. Compare its phase with −ωT for ωT = 0.5, 1, 2 and 3, and at ω_c for the three delays. (d) What delay would bring 180° + phase at crossover to zero?

Answer form: a short proof, a frequency and four phases, a table, a delay.
