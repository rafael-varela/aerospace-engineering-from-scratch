# P02.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a sensor mount by hand)

A 1.2 kg sensor sits on a mount of stiffness 40 kN/m with damping 15 N s/m and is shaken by a 5 N force at f = 25 + s/5 Hz.

(a) Find ω_n, ζ and the frequency ratio r. (b) Find the magnification, the steady amplitude in μm and the phase lag in degrees from (3). (c) Find the undetermined coefficients A and B from (1) and show they are consistent with (b). (d) Give the resonant magnification 1/(2ζ), the peak frequency, the half-power bandwidth in Hz and Q. Confirm everything with `magnification`, `phase_lag`, `undetermined_coefficients`, `half_power_bandwidth` and `quality_factor`.

Answer form: ten numbers with units, one consistency check.

## 2. (Derivation: two routes to the steady state, and the exceptional case)

(a) Substitute x_p = A cos ωt + B sin ωt into m ẍ + c ẋ + k x = F₀ cos ωt and derive the 2×2 system (1); solve it for A and B.
(b) Solve the complex equation m z̈ + c ż + k z = F₀ e^{iωt} with z = Z e^{iωt}, take the real part, and show A = F₀ Re H and B = −F₀ Im H.
(c) Differentiate the denominator of M² with respect to r² to derive r_peak = √(1 − 2ζ²) and M_peak = 1/(2ζ√(1 − ζ²)); state for which ζ no peak exists.
(d) For c = 0 and ω = ω_n, substitute x_p = t(A cos ω_n t + B sin ω_n t) and show A = 0, B = F₀/(2mω_n).
(e) Show that the beats formula (6) tends to (5) as ω → ω_n.

Answer form: five derivations.

## 3. (Unfamiliar scenario: designing an accelerometer)

A piezoresistive accelerometer is a proof mass on a spring with damping; its output is the relative displacement z, which obeys m z̈ + c ż + k z = −m ÿ for base acceleration ÿ. For a sinusoidal base acceleration of amplitude a at frequency ω, the output amplitude is (a/ω_n²) M(r, ζ), so the instrument reads correctly only where M ≈ 1.

(a) Explain in two sentences why an accelerometer needs f_n *above* its measurement band while a vibration isolator needs f_n *below* the disturbance band. (b) The instrument must read within ±5 % up to f_max = 100 + 10s Hz with ζ = 0.7. Find the largest r at which |M − 1| ≤ 0.05 (numerically, from `magnification`) and hence the minimum f_n. (c) Repeat for ζ = 0.5, 0.6, 0.65 and 0.75 and find which ζ allows the lowest f_n for the same ±5 % band. Explain the result from the expansion M⁻² = 1 + (4ζ² − 2)r² + r⁴. (d) Why would a designer not simply choose f_n ten times f_max?

Answer form: an explanation, one r and one f_n, a table of five (ζ, r_max) pairs with a reason, one sentence.

## 4. (Numerical, uses your code: the tolerance band)

The lesson's bracket (f_n = 55 Hz, F₀ = 20 N) must keep its steady amplitude below X_lim = 0.3 + s/50 mm.

(a) With ζ = 0.02, find the band of drive frequencies in which the steady amplitude exceeds X_lim (use `brentq` on `magnification` either side of the peak). (b) Repeat for ζ = 0.05 and ζ = 0.10; for which ζ does the limit cease to be exceeded at any frequency? (c) Compare each band with the half-power bandwidth and explain why the two differ. (d) Starting from rest at the resonant frequency, how long does the amplitude take to reach X_lim? Read it off `full_response` and compare with (F₀|H|)(1 − e^{−ζω_n t}) = X_lim.

Answer form: three bands in Hz, one critical ζ, a comparison, one time in seconds.

## 5. (Judgment: what the sweep and the force model hide)

(a) The engine force on the bracket comes from rotor unbalance, F₀ = m_e e ω². Rewrite the steady amplitude for this force and find its limit as ω → ∞. How does this change the conclusion that "stiffening always helps"?
(b) A test engineer sweeps a shaker from 45 to 65 Hz in 2 s and reports Q = 18 for the bracket. Using part 4 of the lesson, say whether you believe the number, which way it is biased, and what sweep rate you would require.
(c) The model's damping is a single viscous c. Name two physical damping mechanisms in a bolted bracket, say how each depends on amplitude, and state what a shaker test at two force levels would show for each.
(d) Which single measurement would tell you that the bracket has a second mode inside the 50–60 Hz band?

Answer form: a formula with a limit and a conclusion; a judgment with a rate; two mechanisms; one measurement.

## 6. (Hard: a chirp through resonance)

Drive the bracket (ζ = 0.02, F₀ = 20 N) with F(t) = F₀ cos(2π(45t + 10t²/T)), a linear sweep from 45 Hz to 65 Hz over T seconds, using `rk4_forced` from rest with 10 000 steps per second.

(a) Show that the instantaneous frequency of this force is 45 + 20t/T Hz and the sweep rate is 20/T Hz/s. (b) For T = s/5, s, 5s and 20s seconds, record the peak |x| and the instantaneous frequency at which it occurs, and express the peak as a fraction of the steady resonant amplitude F₀/(cω_n). (c) Derive a time-scale argument for the fraction: the amplitude builds as 1 − e^{−t/τ} with τ = 1/(ζω_n), while the drive stays inside the half-power band 2ζf_n for a time Δt = 2ζf_n/(rate). Compare 1 − e^{−Δt/τ} with your measured fractions and comment. (d) Explain the direction of the frequency shift of the apparent peak, and predict what a *downward* sweep would do to it. Verify.

Answer form: a derivation, a four-row table, a comparison, a prediction with a check.
