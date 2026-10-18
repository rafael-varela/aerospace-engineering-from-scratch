# P06.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which results are exact and which are numerical, and label every random time history SYNTHETIC.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: reading an isolator)

An IMU of mass m = 0.30 + s/100 kg sits on four elastomer mounts of 12 kN/m each, with ζ = 0.1 (treat as viscous).

(a) Find fₙ, the damping coefficient c, and the static sag under 1 g and under 5 g. Check fₙ from the 1 g sag with fₙ = (1/2π)√(g/δ).
(b) Find T at r = 1, at r = 2 and at 500 Hz, and the isolation efficiency 1 − T at 500 Hz.
(c) Compare the exact resonant peak with 1/(2ζ). Where is it?

Answer form: three numbers and a check; three transmissibilities and an efficiency; a peak, its frequency ratio and a comparison.

## 2. (Derivation)

(a) Starting from m ẍ = −k(x − y) − c(ẋ − ẏ), derive the transmissibility T(r) and the phase lag. Show that the force transmitted to a rigid base by a mass driven by F₀e^{iωt} has the same ratio |F_T/F₀|.
(b) Prove that T(√2) = 1 for every ζ, and that T < 1 for all r > √2.
(c) Derive the peak location r_p² = (√(1 + 8ζ²) − 1)/(4ζ²), and show T_peak → 1/(2ζ) as ζ → 0.
(d) Derive T for a hysteretic mount with complex stiffness k(1 + iη), and its high-frequency asymptote. Why does viscous damping give −20 dB/decade and hysteretic damping −40 dB/decade?

Answer form: four derivations.

## 3. (Unfamiliar scenario: a camera on a helicopter)

A 2.5 kg camera is mounted on a helicopter whose four-bladed rotor turns at 5.6 rev/s. The dominant airframe vibration is at the blade-passing frequency (4/rev); there is also a 1/rev component from rotor imbalance. Take viscous ζ = 0.1.

(a) What is the highest fₙ that keeps T ≤ 0.15 at 4/rev? What stiffness and 1 g sag does that need?
(b) The camera's stabilisation system tolerates at most 10 % amplification at 1/rev. What is the lowest fₙ that allows? Can one isolator meet both requirements?
(c) If you accept the isolator of (a), what is the amplification at 1/rev and the transmissibility at 8/rev?
(d) Does a hysteretic mount (η = 0.2) change the conclusion? Suggest a different device for the 4/rev problem.

Answer form: a frequency, a stiffness and a sag; a frequency and a verdict; two transmissibilities; a frequency and a suggestion.

## 4. (Numerical, uses your code: your own design)

Repeat the lesson's design with a sway limit of (0.6 + s/100) mm, keeping the other requirements (fₙ ≥ 60 Hz, peak ≤ 4, T(500 Hz) ≤ 0.06, 5 g static).

(a) Search fₙ = 40–120 Hz in 1 Hz steps and ζ = 0.05–0.30 in steps of 0.01 with `choose_isolator`. Report fₙ, ζ, k, c and the four performance numbers. Which requirements bind?
(b) Check the transmitted Grms and the sway RMS with a 32 s SYNTHETIC simulation using your seed (plus 0.5 s discarded). Express each difference in units of the predicted scatter 1/(2√(ζωₙT)).
(c) Compare with the lesson's 5 Hz-grid answer and explain the difference.

Answer form: a design table and a sentence; two normalised deviations; a sentence.

## 5. (Judgment: "add more damping to be safe")

A colleague proposes doubling the damping of the lesson's isolator, from ζ = 0.2 to 0.4, at fₙ = 65 Hz, "to be safe".

(a) Compute the transmitted Grms on the workmanship spectrum, the peak and T(500 Hz) at ζ = 0.1, 0.2 and 0.4. Which requirement does ζ = 0.4 break?
(b) Now take an illustrative high-frequency-dominated input: 0.001 g²/Hz from 20 to 300 Hz, rising to 0.1 g²/Hz at 500 Hz and flat to 2000 Hz. Compute its overall Grms and the transmitted Grms at fₙ = 65 Hz for ζ = 0.05, 0.1, 0.2, 0.4.
(c) When does more damping help, and when does it hurt? Answer in terms of where the input's energy lies relative to fₙ, and give your recommendation to the colleague.

Answer form: a table and a verdict; a Grms and four numbers; a short argument with a recommendation.

## 6. (Numerical and judgment: viscous or hysteretic?)

A low-level sine-sweep test of a mount gives transmissibility readings. Make a SYNTHETIC version: at 25 log-spaced frequencies from 20 to 1000 Hz, T = T_h(f/70 Hz, η = 0.3) × (1 + 0.03 ε), with ε standard normal from `numpy.random.default_rng(s)`.

(a) Fit both the viscous model (fₙ, ζ) and the hysteretic model (fₙ, η) by weighted least squares with `curve_fit` (σ = 0.03 T). Report the parameters with uncertainties and each χ² against its degrees of freedom.
(b) Fit a straight line to log T against log f above 300 Hz. What slope do you get, and what does each model predict?
(c) Which model do the data support, and what would you tell the designer who used ζ = η/2 in a viscous analysis?

Answer form: two fits with χ²; a slope and two predictions; a decision and a sentence.
