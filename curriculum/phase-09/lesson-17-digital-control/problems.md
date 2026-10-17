# P09.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one servo pole, three discretisations)

An elevator servo is modelled as G(s) = 1/(τs + 1) with τ = 0.05 s. The flight computer samples at T = (5 + s) ms.

(a) By hand, give the ZOH equivalent G_d(z) = b/(z − a) with a = e^{−T/τ} and b = 1 − a. (b) Give the Tustin equivalent and its pole (1 − T/(2τ))/(1 + T/(2τ)), and the matched pole-zero equivalent. (c) Check all three with `c2d_zoh`, `tustin_tf` and `matched_tf`, and give the DC gain of each. (d) Repeat (b) for T = 0.15 s. What is the sign of the Tustin pole, and what does a negative real pole do to a step response?

Answer form: three transfer functions with numbers, three DC gains, one pole with a sign and one sentence.

## 2. (Derivation: the four formulas of the lesson)

(a) Derive the augmented-matrix identity exp([[A, B], [0, 0]]T) = [[e^{AT}, ∫₀ᵀ e^{As}ds B], [0, I]] from the power series. (b) For a computational delay 0 ≤ τ ≤ T, derive x_{k+1} = e^{AT}x_k + G₁u_{k−1} + G₀u_k with G₀ = ∫₀^{T−τ} e^{As}ds B and G₁ = e^{A(T−τ)} ∫₀^τ e^{As}ds B, and show G₀ + G₁ = B_d. (c) Show that s = (2/T)(z − 1)/(z + 1) maps the imaginary axis onto the unit circle and the left half plane into the unit disk, and derive the warping law ω_d = (2/T) arctan(ω_a T/2). (d) Show that the zero-order hold has the frequency response H(jω) = T e^{−jωT/2} sinc(ωT/2) (with sinc x = sin x / x) relative to an ideal sampler, so it adds a phase lag of ωT/2.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a notch filter for a bending mode)

A launch vehicle's first bending mode appears in the pitch-rate gyro at f_b = (8 + s/5) Hz. The flight software must remove it with a notch filter N(s) = (s² + 2ζ_z ω_b s + ω_b²)/(s² + 2ζ_p ω_b s + ω_b²), ζ_z = 0.02, ζ_p = 0.5, ω_b = 2πf_b, running at 100 Hz.

(a) Discretise N by plain Tustin, by Tustin prewarped at ω_b, and by matched pole-zero. (b) For each, find the frequency (in Hz) of the deepest attenuation on the unit circle and the attenuation at f_b in dB. (c) Which method would you fly, and why does plain Tustin miss? (d) What sample rate would make plain Tustin's notch fall within 1 % of f_b?

Answer form: three filters, a three-row table, one choice with a reason, one sample rate.

## 4. (Numerical, uses your code: your flight computer)

Your pitch loop runs at f_s = (10 + s) Hz.

(a) Give the spectral radius, the equivalent damping ratio of the least-damped complex pole and the phase margin from `sample_rate_study` and `discrete_loop_margin` for delays τ = 0, T/2 and T. (b) Compare each margin with the rule PM_d ≈ PM_c − ω_c(T/2 + τ) (PM_c = 65.7°, ω_c = 7.12 rad/s). (c) Find, to 1 Hz, the lowest sample rate for which the loop with a full-sample delay keeps PM ≥ 45°. (d) Confirm one case by `hybrid_simulation`: report the largest difference from the difference equation at the samples.

Answer form: a three-row table with the rule column, one sample rate, one difference.

## 5. (Judgment: emulation versus direct digital design)

(a) The lesson designs in continuous time and then samples ("emulation"). Using the margin table, state a rule of thumb for the ratio ω_s/ω_c (ω_s = 2π/T) above which emulation costs less than 10° of phase margin with a full-sample delay. (b) Name two effects in a real flight computer, besides the hold and the computational delay, that add phase lag in the same loop, and estimate one of them (for example a second-order anti-aliasing filter at 40 Hz from `P09.06`). (c) Timing jitter makes τ vary from sample to sample between 0 and T. Which value of τ would you design for, and why? (d) When would you abandon emulation and design the controller directly in the z-domain?

Answer form: one ratio, two effects with one number, one choice with a reason, one paragraph.

## 6. (Hard: the digital gain margin)

A continuous pitch loop of this kind has a large gain margin; sampling and delay make it finite.

(a) For f_s = 20 Hz and f_s = 50 Hz, with τ = 0 and τ = T, scale the gain vector K by k and find (by bisection on the spectral radius of `discrete_pitch_loop(T, tau, k K)`) the critical k at which the loop becomes unstable. (b) At the critical k, give the frequency of the critical pole pair, ln(z)/T, in rad/s. (c) Compare each critical frequency with the Nyquist frequency π/T, and explain why adding the full-sample delay lowers both the critical gain and the critical frequency. (d) Compare with the continuous loop: is there a finite critical k without sampling? (Check the root locus of A − kBK numerically for k up to 100.)

Answer form: a four-entry table of critical gains, four frequencies, a short explanation, a yes/no with evidence.
