# P06.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Label every result derived from a SYNTHETIC record as such.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Material data (typical, not design allowables): 6061-T6 E = 68.9 GPa, ρ = 2700 kg/m³, yield 276 MPa.

## 1. (Warm-up: your own bracket)

A star-tracker baffle of m_tip = (0.4 + s/50) kg sits at the end of a 6061-T6 cantilever bracket L = (80 + s) mm long, 40 mm wide and 5 mm thick.

(a) The bracket's I, k, own mass and the effective mass (2).
(b) The natural frequency with and without the bracket's own mass, and the difference in per cent.
(c) The static tip deflection and root stress under 1 g. How many g of base acceleration at resonance, with ζ = 0.02, would bring the root to yield?

Answer form: four numbers; two frequencies and a percentage; a deflection, a stress and a number of g.

## 2. (Derivation: the pieces of the lesson)

(a) Derive the effective mass m_tip + (33/140) m_beam by Rayleigh's method from the static deflection shape w(x) = w_tip (3x²L − x³)/(2L³).
(b) From m ẍ = −c(ẋ − ẏ) − k(x − y), derive the relative transmissibility (7) and the absolute transmissibility (8). Show that (8) equals 1 at r = √2 for every ζ.
(c) For a base sine at exactly r = 1 switched on at t = 0 with z(0) = ż(0) = 0, show from (10) that for small ζ the envelope of z grows as (A/(2ζω_n²))(1 − e^{−ζω_n t}), and find the time to reach 99 %.
(d) Show that Coulomb friction reduces the amplitude by 4F/k per cycle without changing the period, and derive the equivalent viscous damping ratio (12) by equating the energy dissipated per cycle.

Answer form: four derivations.

## 3. (Unfamiliar scenario: designing a vibration sensor)

The relative motion z of a small mass on a spring inside a case bolted to a vibrating structure is what a vibration sensor measures.

(a) **Accelerometer.** Show that below resonance z ≈ −ÿ/ω_n², so z measures the base acceleration. For the "gain" ω_n²|Z|/|A| = 1/√((1 − r²)² + (2ζr)²) to stay within 5 % of 1, find the largest r for ζ = 0.02, 0.6 and 0.7, and the smallest sensor natural frequency that measures 2 kHz flat to 5 % in each case. What is the phase lag at that r?
(b) **Seismometer.** Show that far above resonance z ≈ −y, so z measures the base displacement. For |Z/Y| (7) within 5 % of 1, find the smallest r for ζ = 0.02 and 0.7.
(c) Which damping would you choose for an accelerometer meant to measure launch transients (with a broad, rapidly changing spectrum)? Why does phase matter?

Answer form: three r values, three frequencies and three phases; two r values; a choice and a reason.

## 4. (Numerical, uses your code: a SYNTHETIC tap test and a test prediction)

Use the lesson's camera bracket.

(a) Generate the SYNTHETIC record `synthetic_tap_test(seed=s, fn=78.27, zeta=0.02)`. Estimate ζ by the log decrement over n = 1, 5, 10, 15 and 20 cycles, and the damped frequency from the mean of the first 10 peak intervals. Which estimate do you trust, and why?
(b) With your best ζ, predict the camera's peak acceleration (g), the bracket's deflection and its root stress for a 2 g resonant input.
(c) Simulate a 1 s resonant dwell from rest with `simulate_base_sine` and compare the largest |z| with (b). What fraction of the steady value has it reached?
(d) Your estimate of ζ is uncertain by the scatter of experiment 1 (about ± 0.0004 at 10 cycles). What range of root stress does that imply?

Answer form: five ζ values, a frequency and a sentence; three numbers; a deflection and a fraction; a stress range.

## 5. (Numerical, harder: friction and viscous damping together)

A joint has both viscous damping ζ = 0.01 and friction F/k = d = 1 μm.

(a) Modify `simulate_coulomb` to include the viscous term and compute the free decay of the lesson's bracket (ω_n = 491.82 rad/s) from x₀ = 100 μm at rest for 40 periods.
(b) Find the positive peaks. Fit them with a straight line (pure friction) and with an exponential (pure viscous). Which fits better? Is either adequate?
(c) Show that the peaks obey x_{n+1} = q x_n − D for constants q and D. Fit q and D by linear regression of x_{n+1} on x_n, recover ζ from q, and show that D = d(1 + √q)². Derive that formula by following one half cycle at a time.
(d) When does the motion stop?

Answer form: a plot; two residuals and a sentence; q, ζ, D and a derivation; a time.

## 6. (Judgment: before the test)

Give a short argued answer to each.

(a) The test lab proposes to notch the input at 78 Hz to keep the camera below 30 g. What base acceleration does that require (with ζ = 0.02), and what is the risk of notching on the basis of an analysis alone?
(b) The tap test was done with a fingertip at about 0.1 g; the test is at 2 g. Give one reason the damping at test level could be lower and one reason it could be higher.
(c) The bracket's factor of 2.2 on yield is for the nominal root stress. Name two things that make the real peak stress higher.
(d) Is the effective-mass model (2) accurate enough here? What would you do to find out (see `P06.15`)?

Answer form: four labelled items.
