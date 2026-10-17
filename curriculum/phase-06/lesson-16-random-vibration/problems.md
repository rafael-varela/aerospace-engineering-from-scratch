# P06.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are approximations and which are statistical estimates (with their scatter). All random time histories are SYNTHETIC; say so in every table and plot.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a qualification-style spectrum)

An illustrative qualification-style spectrum rises at +6 dB/oct from 20 Hz to 50 Hz, is flat at G_p = 0.10 + s/200 g²/Hz from 50 Hz to 800 Hz, and falls at −6 dB/oct to 2000 Hz.

(a) Find G at 20 Hz and at 2000 Hz as multiples of G_p.
(b) Compute the area of each segment by hand with the exact power-law formula, and the overall Grms. Check with `psd_area` and `grms`.
(c) Integrate the same four breakpoints with the trapezoid rule on linear axes (`np.trapezoid(G_bp, f_bp)`). By how much does it overestimate Grms, and why?

Answer form: two multiples; three areas and a Grms; a Grms, a percentage and a sentence.

## 2. (Derivation: Miles' equation and its relatives)

(a) Derive the exact area of a log-log segment G = G₁(f/f₁)ⁿ between f₁ and f₂, including the case n = −1. Which dB/oct slope is n = −1?
(b) With D(r) = (1 − r²)² + 4ζ²r², show that I₂ = ∫₀^∞ r² dr/D equals I₀ = ∫₀^∞ dr/D by the substitution r → 1/r. Then, given I₀ = π/(4ζ), show that ∫₀^∞ |H|² df = (π/2) fₙ Q (1 + 4ζ²) and obtain Miles' equation. State each approximation you make.
(c) ★ Prove I₀ = π/(4ζ) for 0 < ζ < 1, by partial fractions in r² or by contour integration.
(d) Derive the relative-displacement RMS for a white base-acceleration PSD, z_rms = g₀√(G/(64π³ζfₙ³)), and show it equals g₀ Grms_Miles/ωₙ².

Answer form: four derivations.

## 3. (Unfamiliar scenario: a star-tracker bracket on the falling slope)

A 1.5 kg star-tracker head sits on a bracket with its first mode at fₙ = 900 + 10s Hz and ζ = 0.03, on the lesson's workmanship spectrum (NASA-STD-7001C Table 2).

(a) Compute the response Grms by Miles' equation and by `response_grms`. Which is larger, and why, in terms of the shape of the spectrum around fₙ?
(b) What fraction of the full integral comes from frequencies below 0.7fₙ? What does Miles do with that part?
(c) Give the 3σ acceleration and bracket force, and the expected largest peak in a 120 s test.
(d) The tracker's baffle clears its housing by 0.10 mm. Is the 3σ relative displacement a concern?

Answer form: two Grms and a sentence; a fraction and a sentence; two loads and a peak; a displacement and a verdict.

## 4. (Numerical, uses your code: identify the box from a SYNTHETIC test)

(a) With your seed s, synthesize the workmanship input for T = 4, 16 and 64 s (plus 0.5 s to discard), drive the lesson's box (150 Hz, ζ = 0.05) and compare each RMS with the full integral. Express each difference in units of the predicted scatter 1/(2√(ζωₙT)). Are they consistent?
(b) From the 64 s SYNTHETIC input and response, estimate the transfer function with `scipy.signal.csd` and `scipy.signal.welch` (H = P_xy/P_xx, `nperseg=8192`), and fit fₙ and ζ by least squares of |H| on 60–400 Hz with `scipy.optimize.curve_fit`. Compare with the half-power bandwidth estimate of ζ.
(c) Which estimate would you trust on a real shaker test, and what limits the half-power method here?

Answer form: three ratios and three normalised deviations; two fitted values with uncertainties and one bandwidth estimate; a short argument.

## 5. (Judgment: the test came in low)

On the shaker, the box of the lesson responds at 7.2 Grms instead of the predicted 9.66 Grms. The control accelerometer confirms the input spectrum within a control tolerance of ±1.5 dB on the PSD (a typical value; NASA-STD-7001C 4.4.4 sets the tolerances for a real test).

(a) What input shortfall in dB would explain 7.2 Grms on its own? Is it inside the tolerance?
(b) What Q (at fₙ = 150 Hz) would explain it? What fₙ (at Q = 10)?
(c) List three physical causes, rank them, and name the measurement that separates them.
(d) Should the mount design load be lowered to match the test? Argue both ways and decide.

Answer form: a dB value and a verdict; a Q and a frequency; a ranked list with tests; a decision with reasons.

## 6. (Numerical and judgment: how narrow is the band?)

(a) Synthesize 32.5 s with seed s, drive the box, discard 0.5 s, and find every positive local maximum of the response acceleration. Compare the fraction of peaks above 1σ, 2σ and 3σ with the Rayleigh law exp(−k²/2) for a narrow-band Gaussian process.
(b) Compute the irregularity factor α = m₂/√(m₀m₄) of the response PSD (you will need m₄ = ∫f⁴S df) and the rate of maxima √(m₄/m₂). Explain the difference you found in (a).
(c) For a 60 s test, how many cycles exceed 3σ according to the Rayleigh law, and at what level is one exceedance expected? Relate the second number to Davenport's peak factor.

Answer form: a table of three pairs of fractions; two numbers and an explanation; two numbers and a sentence.
