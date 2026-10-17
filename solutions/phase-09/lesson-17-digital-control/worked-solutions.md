# P09.17  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or checked against `scipy.signal.cont2discrete`; the "independent check" line says how it was verified by a different route. Units: SI; frequencies in rad/s unless given in Hz.

## 1. One servo pole

**Assumptions.** τ = 0.05 s, T = (5 + 15) ms = **0.02 s**; unit DC gain.

**Result.** (a) a = e^{−0.4} = 0.67032, b = 0.32968: **G_d,ZOH = 0.32968/(z − 0.67032)**. (b) T/(2τ) = 0.2: **G_d,Tustin = (0.16667z + 0.16667)/(z − 0.66667)**, pole 0.8/1.2 = 0.66667. Matched: pole e^{−0.4} = 0.67032, zero at −1, gain 0.16484: **G_d,matched = (0.16484z + 0.16484)/(z − 0.67032)**. (c) `c2d_zoh`, `tustin_tf`, `matched_tf` reproduce these; DC gains G_d(1) = 1.0000 for all three (ZOH and Tustin by construction for a stable first-order lag, matched by the scaling). (d) T = 0.15 s: T/(2τ) = 1.5, Tustin pole (1 − 1.5)/(1 + 1.5) = **−0.2**, G_d = (0.6z + 0.6)/(z + 0.2). A negative real pole contributes (−0.2)ᵏ: a component that flips sign every sample, an oscillation at the Nyquist frequency the real servo never shows. ZOH at the same T has pole e^{−3} = 0.0498 > 0: no such artefact.

**Independent check.** `cont2discrete(([1], [0.05, 1]), 0.02, method="bilinear")` gives the same Tustin coefficients (test 3 checks this for the lead compensator).

**Common mistakes.** Forgetting the (z + 1) factor in the Tustin numerator; matching the gain at z = −1 instead of z = 1; concluding Tustin is "unstable" at large T (|−0.2| < 1: it is stable, just wrong).

## 2. Derivations

**(a)** M = [[A, B], [0, 0]], M² = [[A², AB], [0, 0]], and by induction Mᵏ = [[Aᵏ, A^{k−1}B], [0, 0]] (k ≥ 1). So e^{MT} = I + Σ_{k≥1} MᵏTᵏ/k! = [[e^{AT}, Σ_{k≥1} A^{k−1}TᵏB/k!], [0, I]], and Σ_{k≥1} A^{k−1}Tᵏ/k! = ∫₀ᵀ Σ_{j≥0} Aʲsʲ/j! ds = ∫₀ᵀ e^{As}ds.

**(b)** On [t_k, t_k + τ) the input is u_{k−1}: x(t_k + τ) = e^{Aτ}x_k + (∫₀^τ e^{As}ds)Bu_{k−1}. On [t_k + τ, t_{k+1}) it is u_k: x_{k+1} = e^{A(T−τ)}x(t_k + τ) + (∫₀^{T−τ} e^{As}ds)Bu_k. Substituting: **x_{k+1} = e^{AT}x_k + G₁u_{k−1} + G₀u_k** with G₁ = e^{A(T−τ)}∫₀^τ e^{As}ds B and G₀ = ∫₀^{T−τ} e^{As}ds B. Then G₁ = ∫₀^τ e^{A(T−τ+σ)}dσ B = ∫_{T−τ}^T e^{As}ds B, and **G₀ + G₁ = ∫₀ᵀ e^{As}ds B = B_d**.

**(c)** With z = e^{jθ}: (z − 1)/(z + 1) = (e^{jθ/2} − e^{−jθ/2})/(e^{jθ/2} + e^{−jθ/2}) = j tan(θ/2). So s = (2j/T)tan(θ/2) is on the imaginary axis: the unit circle maps onto it. Inverting, z = (1 + sT/2)/(1 − sT/2); |z| < 1 ⇔ |s + 2/T| < |s − 2/T| ⇔ Re s < 0. With s = jω_a and θ = ω_dT: **ω_a = (2/T)tan(ω_dT/2), i.e. ω_d = (2/T)arctan(ω_aT/2)**. Prewarping uses c = ω₀/tan(ω₀T/2) in place of 2/T so that ω_d = ω₀ when ω_a = ω₀.

**(d)** The hold turns a unit sample into a pulse of height 1 and width T: H(s) = (1 − e^{−sT})/s. Factor: H = e^{−sT/2}(e^{sT/2} − e^{−sT/2})/s, and at s = jω, (e^{jωT/2} − e^{−jωT/2})/(jω) = 2 sin(ωT/2)/ω = T sinc(ωT/2). So **H(jω) = T e^{−jωT/2} sinc(ωT/2)**: gain T (which cancels the 1/T of the sampler) times a droop, and a pure phase lag of ωT/2.

**Common mistakes.** In (b), forgetting to propagate the τ-interval contribution through e^{A(T−τ)}; in (c), using tan(θ) instead of tan(θ/2).

## 3. Notch filter

**Assumptions.** f_b = 8 + 15/5 = **11 Hz**, ω_b = 69.115 rad/s, f_s = 100 Hz, ζ_z = 0.02, ζ_p = 0.5.

**Result.** (a)

- plain Tustin: N_d = (0.77355z² − 1.20216z + 0.75468)/(z² − 1.20216z + 0.52822);
- prewarped at ω_b: (0.76798z² − 1.16858z + 0.74865)/(z² − 1.16858z + 0.51663);
- matched pole-zero: (0.73224z² − 1.11304z + 0.71228)/(z² − 1.16952z + 0.50100).

(b)

| method | deepest attenuation at | depth (dB) | attenuation at 11 Hz (dB) |
|---|---|---|---|
| continuous | 11.000 Hz | −27.96 | −27.96 |
| plain Tustin | 10.591 Hz | −28.0 | −20.83 |
| Tustin prewarped at ω_b | 11.000 Hz | −28.0 | −27.96 |
| matched pole-zero | 11.000 Hz | −28.0 | −27.96 |

(c) Fly prewarped Tustin (or matched pole-zero; both are exact at the notch here). Plain Tustin moves the notch 3.7 % low, to (2/T)arctan(ω_bT/2)/(2π) = 10.591 Hz, and a lightly damped zero pair is narrow, so the mode at 11 Hz gets 20.8 dB instead of 28 dB: 2.3 times more of the bending signal reaches the controller. (d) Scanning f_s: 150 Hz gives a notch 1.7 % low, **200 Hz** gives 0.98 % low, the first of the scan within 1 %.

**Independent check.** The arctangent law predicts 10.59104 Hz at 100 Hz; the grid minimum is 10.59103 Hz.

**Common mistakes.** Prewarping at the crossover of the control loop instead of at the notch; judging the filter by its depth (all three are −28 dB deep) instead of by where the depth is.

## 4. Your flight computer (s = 15: f_s = 25 Hz, T = 0.04 s)

**Result.** (a)–(b)

| τ | spectral radius | ζ_eq | PM measured (°) | PM rule (°) |
|---|---|---|---|---|
| 0 | 0.9680 | 0.614 | 57.7 | 57.6 |
| T/2 | 0.9678 | 0.495 | 49.6 | 49.4 |
| T | 0.9677 | 0.373 | 41.5 | 41.3 |

The rule is within 0.2° at 25 Hz. (c) Scanning f_s in 1 Hz steps with a full-sample delay: 29 Hz gives 44.76°, **30 Hz gives 45.50°**. (d) `hybrid_simulation(0.04, 0.02, 3.0)` against the difference equation of `discrete_pitch_loop(0.04, 0.02)`: largest difference at the samples **7.5 × 10⁻¹⁵**.

**Common mistakes.** Reading the spectral radius as a margin (it barely changes across the table while the margin loses 16°).

## 5. Judgment

**(a)** With a full-sample delay the loss is ≈ 1.5ω_cT = 540°/(ω_s/ω_c). Ten degrees needs **ω_s/ω_c ≈ 54**. Measured: at 60 Hz (ω_s/ω_c = 52.9) the loss is 10.09°, at 80 Hz (70.6) 7.6°, so about 55 or more. With no computational delay the requirement falls by a factor 3, to about 18.

**(b)** The anti-aliasing filter and the sensor's own bandwidth. A second-order Butterworth anti-aliasing filter at 40 Hz has a phase of −2.30° at 7.12 rad/s. Others: digital pre-filtering of the gyro (`P09.07`), the elevator servo's electronics, bus latency between computers.

**(c)** Design for τ = T. The margin falls monotonically with τ (part 4), so the worst case sets the guarantee; and a scheduler that sometimes finishes early delivers jitter, which a fixed one-sample delay removes (output at the next tick, always).

**(d)** Abandon emulation when the sample rate is only a few times the bandwidth (ω_s/ω_c below about 10, where the rule fails and the margins vanish), when a lightly damped mode lies near the Nyquist frequency, or when a deadbeat or minimum-time discrete response is wanted. Then design on (Φ, Γ) directly (Extension 2).

## 6. Digital gain margin

**Result.** (a)–(b)

| f_s | τ | critical k | critical frequency (rad/s) | Nyquist π/T (rad/s) |
|---|---|---|---|---|
| 20 Hz | 0 | 7.090 | 25.11 | 62.83 |
| 20 Hz | T | 2.247 | 12.32 | 62.83 |
| 50 Hz | 0 | 16.30 | 41.49 | 157.08 |
| 50 Hz | T | 5.457 | 22.29 | 157.08 |

(c) None of the critical frequencies is near Nyquist (they are 40 %, 20 %, 26 % and 14 % of it): the instability is a low-frequency phase-crossover instability, not a sampling artefact at π/T. The delay adds phase lag ωτ, so the loop phase reaches −180° at a lower frequency; the plant gain is larger there (it falls with frequency), so less extra gain is needed to reach |L| = 1. Hence both the critical gain and the frequency fall with τ.

(d) **No finite critical k without sampling**: over k ∈ [0.01, 100], max Re eig(A − kBK) = −0.0236 s⁻¹. As k grows the fast pair moves to −9.25 ± j∞ (−9.24 ± 109.08j at k = 100, −9.25 ± 346.26j at k = 1000): the root-locus asymptotes stay in the left half plane, so the continuous loop has an infinite gain margin. Sampling and delay make it finite (7.09 at 20 Hz even with no delay).

**Independent check.** At each critical k the critical pole has |z| = 1 to the bisection tolerance and Re(ln z)/T ≈ 10⁻¹⁴.

## Check-yourself answers (lesson.md)

1. A_d = e^{AT}, B_d = ∫₀ᵀ e^{As}ds B; exp([[A, B], [0, 0]]T) = [[A_d, B_d], [0, I]].
2. The discrete step response equals the continuous step response at every sample instant. ZOH has it (a step is a held input); Tustin and matched pole-zero do not.
3. All eigenvalues of M strictly inside the unit circle (spectral radius < 1). Under ZOH, z = e^{sT}.
4. The hold's staircase is on average half a sample late; its frequency response is T e^{−jωT/2} sinc(ωT/2), a phase lag of ωT/2.
5. ΔPM ≈ ω_c(T/2 + τ) = 10 × (0.01 + 0.02) = 0.30 rad ≈ 17°.
6. s = (2/T)(z − 1)/(z + 1); ω_d = (2/T)arctan(ω_aT/2). Prewarping makes the response exact at one frequency ω₀; it does not remove the warping at other frequencies, and it does not add the hold's phase lag to a plant model.
7. To z = −1, the highest frequency the sampled system can represent (the Nyquist frequency), so the discrete filter keeps the high-frequency roll-off of the continuous one and stays causal with the same number of zeros as poles.
8. It is set by the slowest real pole (e^{−0.81T} near 1 for every rate), which hardly moves, while the damping and the margin of the fast oscillatory pair fall steadily; it only signals trouble when the loop is already nearly unstable.
