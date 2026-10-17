# P04.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainty on every numerical answer, and say which results are exact, which numerical and which measured. Any dataset produced by `make_synthetic_pendulum_data` is **SYNTHETIC** and must be labelled so wherever it appears.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Use g = 9.80665 m/s² unless a problem says otherwise.

## 1. (Warm-up: plan a pendulum measurement)

A pendulum of length L = 0.5 + s/50 m will be filmed at 30 fps.

(a) Compute T₀, and the exact period and its relative excess over T₀ at 10° and 20° amplitude.
(b) Compute the frame-quantisation standard uncertainty of one timing, and of the period when N = 1, 5, 10, 20 and 50 swings are timed together.
(c) How many swings must be timed together for the frame term per period to be one tenth of the 10° amplitude correction? Is that practical, given that the amplitude decays during the run?

Answer form: three periods, a table of five uncertainties, one integer and two sentences.

## 2. (Derivation: energy, decrement, bandwidth, the large-amplitude pendulum, the string)

(a) Show that E = ½mẋ² + ½kx² is constant for SHM and equals ½kA².
(b) Show that the ratio of successive peaks of an underdamped oscillator is exp(2πζ/√(1 − ζ²)), and invert it to ζ = δ/√(4π² + δ²).
(c) Show that for small ζ the two frequencies at which X = X_max/√2 are ω₀(1 ± ζ) to first order, so the half-power bandwidth is 2ζω₀ = ω₀/Q.
(d) From energy conservation for a pendulum released from rest at θ₀, show that T = 4√(L/g) ∫₀^{θ₀} dθ/√(2(cos θ − cos θ₀)). With the substitution sin(θ/2) = k sin φ, k = sin(θ₀/2), turn it into 4√(L/g) ∫₀^{π/2} dφ/√(1 − k² sin² φ), and expand to first order in k² to get T ≈ T₀(1 + θ₀²/16).
(e) Derive the wave equation for a string under tension T_s with mass per unit length μ, stating where you use the small-slope assumption, and show that c = √(T_s/μ) has units of speed.

Answer form: five derivations.

## 3. (Unfamiliar scenario: a solar-array boom and a reaction wheel)

A deployed solar-array boom has its first bending mode at f₀ = 0.8 Hz with damping ratio ζ = 0.005 + s/3000. A reaction wheel with a small static imbalance shakes the spacecraft at the wheel's spin frequency.

(a) Compute Q, the frequency of peak response and the amplification at that peak relative to a static load of the same size.
(b) Find the band of wheel speeds (in revolutions per minute) in which the boom's response exceeds ten times the static response. Compare its width with the half-power bandwidth.
(c) After an attitude slew, the boom rings freely. How long does the oscillation take to decay to 1 % of its initial amplitude, and how many cycles is that?
(d) What are the amplifications if the wheel runs at 2f₀ and at f₀/2? Recommend an operating rule for the wheel speed.

Answer form: three numbers, a speed band with its width, a time and a cycle count, two numbers and a one-sentence rule.

## 4. (Numerical, uses your code: the wave solver)

A string with tension T_s = 50 + s N, μ = 2.0 g/m and L = 0.65 m is plucked.

(a) Compute c and the first three standing-wave frequencies.
(b) Release a Gaussian pulse exp(−((x − L/2)/0.03 m)²) from rest and run `wave_equation_fd` with 200 cells for one round trip (t = 2L/c) at Courant numbers 1.0, 0.8 and 0.5. For each, report the maximum error against `dalembert` and the relative spread of `wave_energy`. Explain the pattern.
(c) Release mode 3 at Courant 0.5 and run for 1.3 periods with 50, 100, 200 and 400 cells. Tabulate the maximum error against the exact standing wave and report the observed order.
(d) Run the pulse at Courant 1.02 and 1.05 for one round trip and report the maximum |u|. Why does the instability appear even though the physical solution never exceeds 1?

Answer form: one speed and three frequencies, a 3-row table, a 4-row table with an order, two numbers and two sentences.

## 5. (Data reduction, uses your code: a SYNTHETIC pendulum dataset)

Generate `lengths, times = make_synthetic_pendulum_data(seed=s)`. This is **SYNTHETIC** data (generated at g = 9.81 m/s², θ₀ = 0.15 rad, ten swings, five repeats, 30 fps); treat it as you would your own measurement but label it as synthetic.

(a) For each length, reduce to T ± u_T with `period_from_timings` (frame term u = 0.41/30 s per timing) and to g ± u_g with `g_from_pendulum` (u_L = 2 mm).
(b) Apply the amplitude correction and report the corrected g at each length and their weighted mean with its uncertainty.
(c) Fit T = A L^p with `fit_power_law` and report p and g from A, before and after the amplitude correction.
(d) The generating value is 9.81 m/s². Compute E_n for the weighted mean. Without the amplitude correction, what would E_n have been, and what conclusion would you (wrongly) have drawn?

Answer form: a 6-row table, two numbers with uncertainty, four numbers, two E_n values and one sentence.

## 6. (Judgment: oscillator models in the real world)

For each, say which assumption of this lesson's models fails first, what you would observe, and what you would do. Two or three sentences each.

(a) A 2 m pendulum with a 10 g paper-cone bob.
(b) A pendulum hung from a curtain rail.
(c) A flexible antenna boom whose measured log decrement is 0.02 at 1 mm tip amplitude and 0.08 at 20 mm.
(d) Using the FD solver at Courant 0.8 to predict the arrival time of a sharp shock along a cable.
(e) A wing whose damping ratio in flight test falls from 0.03 to 0.005 as speed rises.

Answer form: five short paragraphs.
