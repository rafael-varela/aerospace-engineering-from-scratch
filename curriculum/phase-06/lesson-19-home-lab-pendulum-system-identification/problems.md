# P06.19  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainties on every numerical answer, name the method behind every uncertainty, and label every SYNTHETIC record, table and plot as SYNTHETIC. If you have measured your own pendulum, use your data wherever a problem allows it and say so.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the prediction for your own pendulum)

A pendulum has L = 0.45 + s/30 m (± 3 mm) to the centre of a 100 g rubber ball 40 mm across, in air (ρ = 1.2 kg/m³, C_d = 0.47); g = 9.81 ± 0.01 m/s².

(a) Predict f₀ ± u. Which input dominates the uncertainty?
(b) Predict the frequency at 10°, 20° and 40°. Between which pairs of amplitudes is the difference detectable from 60 zero crossings timed to about 7 × 10⁻⁵ relative?
(c) Predict b₂ and the equivalent damping ratio at 20° from air drag alone, and the amplitude after 120 s if air drag were the only loss.

Answer form: a frequency with uncertainty and a sentence; three frequencies and a verdict; two numbers and an angle.

## 2. (Derivation)

(a) For θ = A cos ωt, compute the energy removed per cycle by viscous, quadratic and Coulomb damping, and derive dA/dt and ζ_eq(A) for each.
(b) Solve dA/dt = −κA^p for p ≠ 1 and for p = 1. Show that A is linear in t for p = 0 and 1/A is linear in t for p = 2.
(c) Solve dA/dt = −αA − βA² (hint: substitute u = 1/A) and obtain the closed form used by `fit_envelope_vq`.
(d) Show that the two-peak logarithmic decrement over n cycles has standard uncertainty √2 s/n when each peak height has relative uncertainty s, and that ζ ≈ δ/(2π) for small δ.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a wing mode in a ground vibration test)

A ground vibration test taps a wing tip and records the tip acceleration at 256 samples per second for 10 s. Make a SYNTHETIC version: y(t) = `damped_cosine(t, 0.8, 0.015, 2π · 4.2/√(1 − 0.015²), 0.3, 0)` in g, plus normal noise of 0.02 g from `numpy.random.default_rng(s)`.

(a) Fit the damped cosine. Report f_d and ζ with covariance uncertainties, and the residual RMS against the noise.
(b) Estimate ζ by the ln(peak) line and by the two-peak decrement over 20 cycles. Compare all three.
(c) Check the covariance with a 100-sample bootstrap (seed s). Is the covariance trustworthy here, and why is this case different from the lesson's pendulum?
(d) How many cycles does the record contain, and how far has the response decayed by its end? Would a longer record improve ζ?

Answer form: two numbers with uncertainties and an RMS; two estimates and a comparison; a bootstrap sd and an argument; two numbers and a sentence.

## 4. (Numerical, uses your code: how long a record?)

Analyse the committed SYNTHETIC record `projects/p06-pendulum-sysid/data/SYNTHETIC_seed619_swing.csv` truncated to its first 30 s, 60 s and the full 120 s.

(a) For each length, report the NLS ζ and ωₙ with covariance and block-bootstrap (blocks of 60 samples, 100 resamples, seed 2) uncertainties, and the residual RMS.
(b) Explain why the covariance uncertainty falls with record length while the block-bootstrap uncertainty does not, using the residual RMS.
(c) What record length would you choose for this pendulum, and what would you change in the *model* rather than the record?

Answer form: a table of three rows; an explanation; a decision.

## 5. (Numerical, uses your code: tell the damping apart)

(a) With your seed s, generate three SYNTHETIC records released at 30° for 60 s with 0.1° noise and no offset: pure Coulomb (b₀ = 0.02 rad/s²), pure viscous (b₁ = 0.02 s⁻¹) and pure quadratic (b₂ = 0.03 rad⁻¹). Fit the envelope exponent p to each and report p ± u.
(b) Generate a mixed record (b₁ = 0.01 s⁻¹, b₂ = 0.02 rad⁻¹) released at 40° for 150 s with 0.2° noise, and fit both `fit_envelope` and `fit_envelope_vq`. Compare b₁ and b₂ with the truth, and say what the fitted p means for a mixture.
(c) Repeat (b) released at 20° for 120 s. How do the uncertainties and the correlation of b₁ and b₂ change, and why?

Answer form: three exponents; two coefficients, an exponent and a sentence; two coefficients with uncertainties, a correlation and an explanation.

## 6. (Judgment: which uncertainty do you report?)

For the committed SYNTHETIC record, the NLS fit gives ζ = 9.21 × 10⁻⁴ with a covariance uncertainty of 0.10 × 10⁻⁴, a single-residual bootstrap of 0.09 × 10⁻⁴, a block bootstrap of 0.48 × 10⁻⁴, and six repeat runs with fresh noise that scatter by 0.05 × 10⁻⁴. The residual RMS is 0.80° against a known reading noise of 0.30°.

(a) Which uncertainty do you report, and why are the other three smaller?
(b) Fit the damped cosine separately in four consecutive 30 s windows. Report each window's mean amplitude, ζ, residual RMS and amplitude-corrected f₀ (ωₙ divided by the mean ω(A)/ω₀ of its peaks, over 2π), with covariance uncertainties. Is the model adequate in each window? Does ζ follow the amplitude as the generator's ζ_eq(A) says?
(c) Combine the four f₀ values. Is their scatter consistent with their covariance uncertainties? What do you conclude, and what single number for f₀ would you report?

Answer form: a choice with reasons; a table and two verdicts; a weighted mean, a χ², a conclusion and a number with uncertainty.
