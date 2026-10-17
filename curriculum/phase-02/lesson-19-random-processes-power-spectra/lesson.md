# P02.19  Random Processes & Power Spectra

> A gust has no formula, but its statistics do: the autocorrelation says how long it remembers, the spectrum says where its energy sits, and a linear system only needs the spectrum.

| | |
|---|---|
| **Learning outcomes** | 1. Derive, from the autocorrelation function, the one-sided power spectral density of a stationary random process and prove the normalisation ∫Φ dω = σ² for the Dryden gust forms. 2. Implement and test, from scratch, a shaping filter that turns seeded white noise into a Dryden gust record, the biased autocorrelation estimator, the periodogram and Welch's averaged spectrum, against `scipy.signal` as the oracle. 3. Measure the scatter of a variance estimate against record length and show it follows √(5L/(4VT)). 4. Decide, from the RMS load factor and the exceedance fraction of a first-order heave model, whether a light aircraft's ride in moderate low-altitude turbulence is acceptable and what speed change would make it so. |
| **Prerequisites** | `P02.16` (Gaussian random variables, variance, the central limit theorem), `P02.10` (Fourier transform, FFT, sampling and aliasing), `P02.08` (transfer functions, |H(iω)|), `P02.05` (matrix exponential, zero-order-hold discretisation), `P02.03` (first-order lag response), `P01.04` (improper integrals), `P00.03`, `P00.07`. Skill check: you can evaluate ∫₀^∞ dx/(1 + x²) and ∫₀^∞ dx/(1 + x²)², and you know what |H(iω)|² does to a sinusoid's amplitude. |
| **Engineering question** | A 1000 kg light aircraft flies at 60 m/s through low-altitude turbulence with a vertical-gust RMS of 1.5 m/s and a scale length of 300 m: what is the RMS load factor its occupants feel, what fraction of the time does the bump exceed 0.3 g, how often does the plunge velocity reverse, and how long a record is needed to measure the gust variance to 2 %? |
| **Assumptions and limits** | The gust velocity is a stationary, zero-mean, Gaussian random process, frozen in space and swept past the aircraft at constant airspeed V (Taylor's hypothesis), uniform across the span, with the Dryden rational spectrum (its ω⁻² tail is a modelling convenience, not measured physics; real turbulence falls as ω^(−5/3) in the inertial range, which the von Kármán form captures and the Dryden form does not). The aircraft is a plunge-only point mass with quasi-steady lift: no pitching, no unsteady aerodynamics, no structural modes. Estimates assume ergodicity (one long record stands for the ensemble) and records much longer than L/V. The model says nothing about discrete gusts, wind shear, or the non-Gaussian tails of real patches of turbulence. |
| **Analytical work** | Stationarity and the autocorrelation; the Wiener–Khinchin pair and the one-sided convention; normalisation of the Dryden forms; the shaping-filter factorisation |H|² = πΦ; sampled white noise of variance 1/Δt; the heave-mode transfer function and the closed-form response variances by partial fractions; the sharp-edged-gust limit; the variance of a variance estimate; Rice's level-crossing formula and why the Dryden gust has no crossing rate. |
| **Computational deliverable** | `starter/main.py`: `dryden_psd`, `dryden_psd_hz`, `dryden_autocorrelation`, `psd_from_autocorrelation`, `integrated_variance`, `shaping_filter`, `discretize_zoh`, `simulate_lti`, `white_noise`, `generate_gust`, `autocorrelation_estimate`, `periodogram`, `welch_psd`, `band_average`, `heave_time_constant`, `heave_response`, `heave_response_variances`, `response_variance_numeric`, `simulate_heave`, `gaussian_exceedance_fraction`, `rice_upcrossing_rate`, `count_upcrossings`, `variance_estimate_relative_sd`, `record_length_for_variance`. Tests prove: both Dryden spectra integrate to σ² and have the right limits; the autocorrelation transforms to the spectrum; the filter gain squared is πΦ; the discretisation and simulation match `scipy.signal.cont2discrete` and `dlsim`; the seeded gust record has the model variance, mean and Gaussian tails within derived statistical tolerances; the estimators match `scipy.signal.welch`, `periodogram` and `correlate` to 10⁻¹⁰; the Welch estimate reproduces the spectrum band by band and the ω⁻² slope; the periodogram's scatter equals its mean; the closed-form response variances match quadrature and their limits; the simulated heave response, crossing counts and exceedance fraction match the formulas; and the variance scatter follows the 1/√T law. |
| **Practical activity** | **computational experiment**: generate seeded synthetic gust records of 100, 400 and 1600 s, measure the scatter of the variance estimate over 24 seeds against the derived law, and compare periodogram and Welch estimates of the same spectrum. All data in this lesson are synthetic. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `gust_spectrum.png`, the hand derivation of the normalisation and of the response variances, the scatter table, and a short note giving the RMS load factor, the exceedance fraction and the record-length recommendation with the assumption you trust least. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the normalisation proof and the σ_v² derivation are complete; all public tests pass; the hand values for the lesson aircraft agree with the code to three figures; the note states a decision with a number and a range of validity. |
| **References** | `@mit-ocw-6-011-signals-systems-and-inference` Lecture 15 "Normal Equations, Random Processes", Lecture 16 "Wide-sense Stationary Processes, LTI Filtering of WSS Processes", Lecture 17 "LTI Filtering of WSS Processes", Lecture 18 "Power Spectral Density", Lecture 19 "Einstein-Wiener-Khinchin Theorem, PSD Applications, Modeling Filters" (required-free). `@mit-ocw-6-011-signals-systems-and-inference` Lecture 12 "Probabilistic Models, Random Variables", Lecture 13 "Vector Picture for First- and Second-Order Statistics, MMSE and LMMSE Estimation" (supplementary-free, second-order statistics). `@mit-ocw-6-003-signals-and-systems` Lecture 16 "Fourier Transform", Lecture 21 "Sampling" (prerequisite refresher). The Dryden turbulence forms themselves are not in these sources; part 2 defines and verifies them. |
| **Extensions** | ★ Implement the von Kármán spectrum (ω^(−5/3) tail), fit a rational approximation to it, and compare σ_n with the Dryden value. ★ Add the short-period pitch mode (P02.04) to the heave model and show how pitch damping changes σ_n. ★ Generate the u, v and w components together and show they are uncorrelated at every lag. |
| **Time** | unmeasured |

## 1. The engineering problem

Every aircraft spends part of its life in turbulent air. The gusts are not signals you can write down: two flights through the same valley meet different gusts, and a third flight meets different ones again. Yet the loads they cause are the same in a statistical sense, and it is the statistics, not any particular gust, that the structure and the autopilot are designed against. The vertical gust velocity w_g(t) seen by an aircraft is modelled as a **random process**: a family of possible time histories with a probability law. What can be specified is its mean (zero), its variance σ_w², and how its values at two instants are related, the **autocorrelation function** R(τ). The Fourier transform of R is the **power spectral density** Φ(ω), which says how much of the variance lives at each frequency. For low-altitude turbulence the standard engineering choice is the **Dryden** spectrum: a rational function of frequency with two parameters, the RMS σ_w and the scale length L, that appears in the flying-qualities specifications used for aircraft certification and simulator work.

The payoff of describing a gust by its spectrum is this: when a random process passes through a linear time-invariant system with frequency response H(iω), the output spectrum is |H(iω)|² times the input spectrum. The aircraft's response to turbulence then reduces to one integral. Here the "system" is the simplest aircraft model there is, a point mass whose lift responds instantly to the relative angle of attack, so that its vertical velocity follows the gust through a first-order lag with time constant τ = 2m/(ρVSC_Lα). For the lesson aircraft (m = 1000 kg, S = 16 m², C_Lα = 5 rad⁻¹, V = 60 m/s, ρ = 1.19 kg/m³ at 300 m) τ = 0.350 s. The turbulence is "moderate": σ_w = 1.5 m/s and L = 300 m.

The questions are four. What is the RMS load factor the occupants feel? What fraction of the time does the bump exceed 0.3 g? How often does the aircraft's plunge velocity reverse (a crossing-rate question, which turns out to expose a flaw in the gust model)? And, if we wanted to measure σ_w from a flight record, how long must the record be for 2 % accuracy? Success criteria, written first: the from-scratch spectral estimators agree with `scipy.signal` to 10⁻¹⁰ on the same record; the seeded synthetic record reproduces σ_w² within the scatter we derive (5.6 % for 2000 s); the Welch estimate matches the Dryden form band by band within that scatter; the closed-form response variances agree with numerical quadrature to 10⁻⁸ and with the simulation within its own scatter; and the counted zero crossings of the plunge velocity agree with Rice's formula within 25 %.

## 2. Understand and derive

### 2.1 Stationary processes and the autocorrelation

A random process x(t) is a rule that assigns a random variable to every instant t. It is **wide-sense stationary** (WSS) when its mean is constant and the correlation between x(t) and x(t + τ) depends only on the lag τ (`@mit-ocw-6-011-signals-systems-and-inference` Lecture 15, 16):

E[x(t)] = μ,   **R(τ) = E[(x(t) − μ)(x(t + τ) − μ)]**.   (1)

(With this definition R is the autocovariance; for the zero-mean gusts of this lesson the two coincide and we call it the autocorrelation.) Three facts follow at once: R(0) = σ² is the variance; R is even, R(−τ) = R(τ), because the pair (t, t + τ) is the pair (t + τ, t) read backwards; and |R(τ)| ≤ R(0), by the Cauchy–Schwarz inequality. A process is **ergodic** when one long record has the same statistics as the ensemble, so that time averages can replace expectations. That is an assumption, not a theorem, and it is what lets a single flight record estimate R: R̂(k Δt) = (1/N) Σₙ x[n] x[n + k].

The lag at which R has fallen to a small fraction of σ² is the **correlation time**: the memory of the process. For turbulence it is set by L/V, the time to fly through one scale length. Units: R is in (m/s)² for a velocity, and it is dimensionally the same as a variance.

### 2.2 The power spectral density and the Wiener–Khinchin theorem

The **two-sided power spectral density** is the Fourier transform of the autocorrelation (`@mit-ocw-6-011-signals-systems-and-inference` Lecture 18, 19; the Einstein–Wiener–Khinchin theorem):

S(ω) = ∫_{−∞}^{∞} R(τ) e^{−iωτ} dτ,   R(τ) = (1/2π) ∫_{−∞}^{∞} S(ω) e^{iωτ} dω.   (2)

Because R is real and even, S is real, even and (this is the theorem's content) non-negative. Setting τ = 0 in the inverse transform gives σ² = (1/2π)∫S dω: the spectrum tells how the variance is spread over frequency, which is why it is called a *power* spectrum. Engineers fold the negative frequencies onto the positive ones and define the **one-sided** PSD. This lesson's convention, used by every function in the code, is

**Φ(ω) = S(ω)/π for ω ≥ 0,   so that   σ² = ∫₀^∞ Φ(ω) dω**,   (3)

and in hertz, G(f) = 2πΦ(2πf) with σ² = ∫₀^∞ G(f) df. Units: Φ is (m/s)² per rad/s, G is (m/s)²/Hz. Conventions differ between books and between software (per rad/s or per Hz, one- or two-sided, a factor 2π placed in the forward or the inverse transform): the only safe habit is to integrate your spectrum and check that you get the variance. Every spectrum in this lesson is checked that way, and `scipy.signal.welch` returns G(f), per hertz and one-sided, which is what `dryden_psd_hz` is for.

### 2.3 The Dryden gust model

The Dryden model describes each gust component by a scale length L and an RMS σ. With V the airspeed, b = V/L is the **break frequency** in rad/s and x = ω/b = Lω/V is a dimensionless frequency. The spectra are

**Φ_u(ω) = σ_u² (2L_u/(πV)) · 1/(1 + x²)**   (longitudinal),   (4)

**Φ_w(ω) = σ_w² (L_w/(πV)) · (1 + 3x²)/(1 + x²)²**   (vertical; the lateral component v has the same form).   (5)

These are the standard Dryden forms as used in flying-qualities specifications, written here in the one-sided convention (3). Documents differ in where the factor 2 sits (some define the vertical scale length as half of ours and write 2L_w in (5)) and in whether the spectrum is one- or two-sided, so we define everything ourselves and prove the normalisation rather than quote it. Substituting ω = bx, dω = b dx:

∫₀^∞ Φ_w dω = (σ_w²/π) ∫₀^∞ (1 + 3x²)/(1 + x²)² dx = (σ_w²/π)(π/4 + 3π/4) = σ_w²,

using ∫₀^∞ dx/(1 + x²)² = π/4 and ∫₀^∞ x² dx/(1 + x²)² = π/4 (P01.04: substitute x = tan θ). For (4), ∫₀^∞ dx/(1 + x²) = π/2 gives (σ_u²/π)·2·(π/2) = σ_u². Both check, and `test_dryden_psd_integrates_to_the_variance_for_both_components` repeats the check numerically with `quad`.

What the forms say. At zero frequency Φ_w(0) = σ_w² L/(πV): for the lesson values 2.25 × 300/(π × 60) = 3.581 (m/s)² per rad/s, and Φ_u(0) is twice that. The u spectrum falls monotonically; the w spectrum rises slightly, to 9/8 of its zero-frequency value at x = 1/√3 (ω = 0.115 rad/s), before falling. Far above the break frequency both fall as ω⁻²: Φ_w → 3σ_w²V/(πLω²) and Φ_u → 2σ_u²V/(πLω²), so at high frequency the vertical gust carries 1.5 times the energy of the longitudinal one. The scale length is 300 m and the aircraft flies it in 5 s, so b = 0.2 rad/s = 0.032 Hz: almost all the gust energy is below 0.1 Hz, far slower than anything a pilot would call a bump. The bumps come from the tail, and the aircraft's response (2.6) decides how much of the tail matters.

The autocorrelations are the inverse transforms. With a = V/L,

**R_u(τ) = σ_u² e^{−a|τ|}**,   **R_w(τ) = σ_w² (1 − a|τ|/2) e^{−a|τ|}**.   (6)

Check the second by transforming forward: the transform of e^{−a|τ|} is 2a/(a² + ω²) and that of |τ|e^{−a|τ|} is 2(a² − ω²)/(a² + ω²)², so S_w = σ_w²[2a/(a² + ω²) − a(a² − ω²)/(a² + ω²)²] = σ_w² a(a² + 3ω²)/(a² + ω²)², and S_w/π is (5). The vertical correlation crosses zero at τ = 2L/V = 10 s and has a shallow negative lobe after it: a vertical gust is, on average, followed by a weaker gust the other way, which is what a divergence-free velocity field must do. The **integral time scale** ∫₀^∞ R dτ/σ² is L/V = 5 s for u and L/(2V) = 2.5 s for w; the test checks both.

```
   R_w(τ)/σ²                               Φ_w(ω)  (log-log)
   1 ┤╲                                    ──────╮
     │ ╲                                          ╲   slope −2
     │  ╲                                          ╲
   0 ┼────╲────────────────── τ                     ╲
     │     ╲____------                               ╲
     │      2L/V  (negative lobe)              b=V/L  ω
```

### 2.4 Where the model stops being valid

Taylor's hypothesis (frozen turbulence swept past at V) holds when V is large compared with the gust velocities, which it is here (60 against 1.5 m/s). The Dryden ω⁻² tail is the price of a rational spectrum: real turbulence in the inertial range falls as ω^(−5/3) (the von Kármán form has it), so Dryden overstates high-frequency energy by a factor that grows as ω^(1/3). Anything sensitive to the tail is suspect: the gust's own derivative has infinite variance in the Dryden model, since ∫ω²Φ_w dω diverges (the integrand tends to 3σ_w²V/(πL), a constant), and we will see in 2.7 that this makes a crossing rate for the gust meaningless. Stationarity fails as terrain, altitude and time of day change; a 33-minute record over varied ground is not one process but several. Gaussianity fails in the tails: real turbulence has patches, and a Gaussian model underestimates rare large gusts. The specifications handle the last point by making σ_w itself random (a probability of encountering each intensity), which is a Bayesian idea you will meet in P02.20.

### 2.5 From white noise to a gust: the shaping filter

**White noise** is the process with R(τ) = S₀ δ(τ): no memory at all, and a flat two-sided spectrum S(ω) = S₀. It has infinite variance, so it is an idealisation, but it is the right input for generating anything else: passing white noise of two-sided PSD S₀ through a stable filter H(s) gives an output with S_out(ω) = |H(iω)|² S₀ (`@mit-ocw-6-011-signals-systems-and-inference` Lecture 16, 17). To generate a process with one-sided PSD Φ we therefore need |H(iω)|² = πΦ(ω) with S₀ = 1, by (3). Because the Dryden spectra are rational in ω², they factor exactly:

**H_u(s) = σ_u √(2L/V) · b/(s + b)**,   **H_w(s) = σ_w √(L/V) · b(√3 s + b)/(s + b)²**,   b = V/L.   (7)

Check H_w: |b(√3 iω + b)|²/|(iω + b)²|² = b²(b² + 3ω²)/(b² + ω²)² = (1 + 3x²)/(1 + x²)², so |H_w|² = σ_w²(L/V)(1 + 3x²)/(1 + x²)² = πΦ_w. The poles are at s = −b (stable) and the zero at s = −b/√3; choosing the stable, minimum-phase factor is what "modelling filter" means in the Lecture 19 title. In state-space form, H_w is the controllable canonical pair

A = [[0, 1], [−b², −2b]],   B = [0, 1]ᵀ,   C = σ_w √(L/V) · b · [b, √3],   D = 0,   (8)

and H_u is the scalar A = −b, B = b, C = σ_u √(2L/V), D = 0. `test_shaping_filter_gain_squared_is_pi_times_the_psd` evaluates C(iωI − A)⁻¹B + D and compares its squared modulus with πΦ at six frequencies.

**Sampled white noise.** A computer cannot make a δ-correlated signal; it makes a sequence e[n] of independent Gaussians with variance σ_e² at spacing Δt. Its autocorrelation is σ_e² at lag 0 and zero elsewhere, which stands in for a continuous R(τ) = σ_e²Δt δ(τ): the area under the discrete spike. For S₀ = 1 we need **σ_e² = 1/Δt** (`white_noise`). The smaller the step, the larger the samples: that is not a bug, it is what a flat spectrum up to the Nyquist frequency π/Δt costs. The filter is then stepped with the **zero-order hold** discretisation of P02.05 (the input held constant over each step), x[n+1] = A_d x[n] + B_d e[n], with [A_d, B_d] read from the matrix exponential of [[A, B], [0, 0]]Δt. The output at the sample instants then has, to first order in bΔt = 0.004, the stationary variance of the continuous process; the exact sampled statistics (a discrete Lyapunov equation) are within 3 × 10⁻⁶ of σ_w² at Δt = 0.02 s.

### 2.6 The aircraft as a filter: the heave mode

Let v be the aircraft's vertical velocity (up positive) and w_g the vertical gust. The relative vertical velocity of the air is w_g − v, the angle-of-attack increment is (w_g − v)/V, and the lift increment with quasi-steady aerodynamics is ½ρV²S C_Lα (w_g − v)/V. Newton's second law for a point mass that may plunge but not pitch:

m dv/dt = ½ρVSC_Lα (w_g − v)   ⟹   **τ dv/dt + v = w_g,   τ = 2m/(ρVSC_Lα)**.   (9)

Units: kg/(kg m⁻³ · m s⁻¹ · m² · 1) = s. For the lesson aircraft τ = 2000/(1.19 × 60 × 16 × 5) = 0.3501 s. The transfer functions from gust to plunge velocity and to vertical acceleration a = dv/dt are

**H_v(iω) = 1/(1 + iωτ)**,   **H_a(iω) = iω/(1 + iωτ)**,   (10)

a first-order lag (P02.03) and its derivative. At low frequency the aircraft rides the gust, v → w_g, and feels no acceleration; at high frequency it cannot follow and |H_a| → 1/τ. The load-factor increment is Δn = a/g, so Δn's spectrum is |H_a|²Φ_w/g². A **sharp-edged gust**, a step of w_g, gives a(0⁺) = w_g/τ: for a one-σ step, 1.5/0.350 = 4.28 m/s² = 0.437 g. This is the classical quasi-static gust load; the spectral answer is smaller because most of the gust energy is at frequencies the aircraft follows.

**Closed-form response variances.** With c = 1/τ, |H_v|² = c²/(ω² + c²) and |H_a|² = c²ω²/(ω² + c²). The output variance is ∫₀^∞ |H|²Φ_w dω, a rational integrand in y = ω². Partial fractions in y, A/(y + c²) + B/(y + b²) + C/(y + b²)², with ∫₀^∞ dω/(ω² + c²) = π/(2c) and ∫₀^∞ dω/(ω² + b²)² = π/(4b³), give after simplification

**σ_v² = σ_w² · c(b + 2c)/(2(b + c)²)**,   **σ_a² = σ_w² · bc²(2b + 3c)/(2(b + c)²)**.   (11)

Limits, which are the checks: c → ∞ (an aircraft that follows instantly) gives σ_v² → σ_w² and σ_a² → ∞, the infinite derivative variance of 2.4; c → 0 (an aircraft too heavy to respond) gives both → 0. For the lesson values, b = 0.200 rad/s and c = 2.857 rad/s: σ_v = 1.426 m/s (95 % of the gust RMS: the aircraft follows almost everything) and σ_a = 1.328 m/s², an **RMS load factor of 0.135 g**, less than a third of the sharp-edged value. `response_variance_numeric` repeats the integral with `quad`; the two agree to 10⁻⁸. Note that with r = c/b = L/(Vτ) ≫ 1 (here 14.3), σ_a² ≈ (3/2)σ_w² bc ∝ V²: the gust load grows linearly with airspeed, as the classical gust-load formula also says.

**Gaussian tails.** A linear filter keeps a Gaussian process Gaussian, so Δn is Gaussian with SD σ_n = 0.135 g, and the fraction of time with |Δn| > n₀ is erfc(n₀/(σ_n√2)): 14.0 % beyond 0.2 g, 2.67 % beyond 0.3 g, 2.2 × 10⁻⁴ beyond 0.5 g (0.8 s per hour).

### 2.7 How often: Rice's crossing rate, and a model that cannot give one

For a stationary Gaussian process x with derivative variance σ_ẋ², the expected number of upward crossings of the level x₀ per unit time is Rice's formula,

**ν(x₀) = (σ_ẋ/(2πσ_x)) exp(−x₀²/(2σ_x²))**.   (12)

The derivation counts the probability that x lies just below x₀ with a positive slope large enough to cross within dt; the Gaussian joint density of x and ẋ (independent at the same instant, because R′(0) = 0) does the rest. The zero-crossing rate ν₀ = σ_ẋ/(2πσ_x) is the "apparent frequency" of the process. For the plunge velocity v, σ_v̇ = σ_a and (12) gives ν₀ = 1.328/(2π × 1.426) = 0.148 s⁻¹: the aircraft's vertical velocity reverses sign 533 times an hour, and crosses +1σ_v upward e^(−1/2) = 0.61 times as often. For the gust w_g itself, or for the acceleration a, the formula needs σ_ẇ or σ_ȧ, and both are infinite in the Dryden model (the spectra have ω⁻² and ω⁰ tails). The model's answer "infinitely many crossings per second" is the model telling you its tail is unphysical. Real aircraft do have finite load-factor crossing rates, because unsteady lift build-up and structural filtering roll the response off; that is why continuous-turbulence load analyses add those dynamics before applying (12), and why this lesson applies it to v only.

### 2.8 How long a record: the variance of a variance estimate

The sample variance σ̂² = (1/T)∫₀^T x² dt of a zero-mean Gaussian process is itself random. For T much longer than the correlation time, Var(σ̂²) = (2/T) ∫_{−∞}^{∞} R²(τ) dτ (the fourth moment of a Gaussian pair is 2R² plus the square of the variance). For R_w, ∫₀^∞ (1 − aτ/2)² e^{−2aτ} dτ = 1/(2a) − 1/(4a) + 1/(16a) = 5/(16a), so

**Var(σ̂_w²) = 5σ_w⁴ L/(4VT),   relative SD = √(5L/(4VT))**.   (13)

For L = 300 m and V = 60 m/s a 2000 s record gives 5.6 %, and 2 % needs 260 minutes of stationary turbulence: longer than any patch of real turbulence stays stationary, which is why flight-test spectra are quoted with error bars that are never small. Equation (13) also sets the test tolerances below: a tolerance of four standard deviations is 22 % on the variance of the 2000 s record.

The same reasoning explains the two spectral estimators of part 3. The **periodogram** |X_k|²Δt/N has, at every bin, a standard deviation equal to its mean, however long the record: it is an inconsistent estimator, because a longer record gives finer frequency resolution instead of less scatter. **Welch's method** cuts the record into K windowed segments, averages their periodograms, and trades resolution for a relative SD of about 1/√K.

## 3. Build the computational model

Open `starter/main.py`. The dictionaries `AIRCRAFT` and `TURBULENCE` hold the lesson values; every function takes its parameters explicitly so the tests can vary them.

### Step 1: the Dryden forms and their transform pair

```python
def dryden_psd(omega, sigma, L, V, component="w"):
    x2 = (L * np.asarray(omega, float) / V) ** 2
    if component == "u":
        return sigma**2 * (2.0 * L / (np.pi * V)) / (1.0 + x2)                     # (4)
    return sigma**2 * (L / (np.pi * V)) * (1.0 + 3.0 * x2) / (1.0 + x2) ** 2         # (5)

def dryden_autocorrelation(tau, sigma, L, V, component="w"):
    s = V * np.abs(np.asarray(tau, float)) / L
    return sigma**2 * np.exp(-s) if component == "u" else sigma**2 * (1.0 - 0.5 * s) * np.exp(-s)   # (6)
```

Validate the inputs (positive σ, L, V; a known component) and raise `ValueError` otherwise. `dryden_psd_hz` is 2πΦ(2πf). `psd_from_autocorrelation(omega, R, tau_max)` evaluates (2/π)∫₀^{τ_max} R(τ) cos(ωτ) dτ with `quad(..., weight="cos", wvar=omega)`, the numerical Wiener–Khinchin check; `integrated_variance` is `quad(psd, 0, np.inf)`.

### Step 2: the shaping filter, its discretisation, and the simulation loop

```python
def shaping_filter(sigma, L, V, component="w"):
    b = V / L
    if component == "u":
        return np.array([[-b]]), np.array([[b]]), np.array([[sigma * np.sqrt(2 * L / V)]]), np.array([[0.0]])
    K = sigma * np.sqrt(L / V) * b
    return (np.array([[0.0, 1.0], [-b * b, -2 * b]]), np.array([[0.0], [1.0]]),
            np.array([[K * b, K * np.sqrt(3.0)]]), np.array([[0.0]]))                   # (8)

def discretize_zoh(A, B, dt):
    n, m = A.shape[0], B.shape[1]
    M = np.zeros((n + m, n + m)); M[:n, :n] = A * dt; M[:n, n:] = B * dt
    E = expm(M)
    return E[:n, :n], E[:n, n:]

def simulate_lti(Ad, Bd, C, D, u, x0=None):
    x = np.zeros(Ad.shape[0]) if x0 is None else np.array(x0, float)
    y = np.empty(len(u))
    for n in range(len(u)):
        y[n] = C @ x + d * u[n]
        x = Ad @ x + Bd * u[n]
    return y
```

`white_noise(n, dt, rng)` draws `rng.normal(0, 1/sqrt(dt), n)`; `generate_gust(duration, dt, sigma, L, V, component, seed)` chains the three, simulating an extra 10 L/V seconds first and discarding them so the record starts in the stationary state. The test compares `discretize_zoh` with `scipy.signal.cont2discrete(method="zoh")` and `simulate_lti` with `scipy.signal.dlsim`, and checks that a unit step into the first-order lag rises as 1 − e^(−t/τ) *exactly* at the sample instants, which is what zero-order hold guarantees.

### Step 3: the estimators

```python
def autocorrelation_estimate(x, max_lag):
    x = x - x.mean(); N = x.size
    return np.array([np.dot(x[:N - k], x[k:]) / N for k in range(max_lag + 1)])

def periodogram(x, dt):
    N = x.size
    P = np.abs(np.fft.rfft(x)) ** 2 * dt / N
    P[1:] *= 2.0
    if N % 2 == 0: P[-1] /= 2.0
    return np.fft.rfftfreq(N, dt), P

def welch_psd(x, dt, nperseg, overlap=0.5):
    step = nperseg - int(round(overlap * nperseg))
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(nperseg) / nperseg)   # periodic Hann
    acc, K = np.zeros(nperseg // 2 + 1), 0
    for s in range(0, x.size - nperseg + 1, step):
        seg = x[s:s + nperseg]; seg = (seg - seg.mean()) * w
        acc += np.abs(np.fft.rfft(seg)) ** 2; K += 1
    P = acc * dt / np.sum(w**2) / K
    P[1:] *= 2.0
    if nperseg % 2 == 0: P[-1] /= 2.0
    return np.fft.rfftfreq(nperseg, dt), P
```

Three details decide whether this matches `scipy.signal.welch` to 10⁻¹⁰ or only to 10 %: the window is the *periodic* Hann (denominator `nperseg`, not `nperseg − 1`); the scale is Δt/Σw² (so that a white input of variance 1/Δt gives P = 2 at every interior bin: one-sided, per hertz, with S₀ = 1); and each segment has its own mean removed before windowing, which is SciPy's default `detrend="constant"`. The 1/N in the autocorrelation (rather than 1/(N − k)) keeps the estimate positive semi-definite at the price of a bias of order k/N. `band_average(f, P, f_lo, f_hi)` averages bins in a band; it is how the tests compare a noisy estimate with a smooth model.

### Step 4: the aircraft

`heave_time_constant` is (9). `heave_response(omega, tau)` returns (H_v, H_a) from (10). `heave_response_variances(sigma_w, L, V, tau)` is (11) and `response_variance_numeric` its `quad` check. `simulate_heave(t, w, tau)` steps (9) with the gust record as input, but **not** with a zero-order hold: a gust held constant between samples jumps at every sample, and because the Dryden gust has a 1/ω² tail those jumps scale as √Δt rather than Δt. The sampled acceleration then carries the jumps and σ_a comes out 5 % high at Δt = 0.02 s (try it). Taking the gust as *linear* between samples (a first-order hold) and solving (9) exactly over one step gives, with h = Δt/τ,

v[n+1] = e^{−h} v[n] + (1 − e^{−h}) w[n] + (1 − (1 − e^{−h})/h)(w[n+1] − w[n]),   a[n] = (w[n] − v[n])/τ,

which is unbiased to within the 0.1 % of energy above the Nyquist frequency. Then `gaussian_exceedance_fraction` is erfc(x₀/(σ√2)), `rice_upcrossing_rate` is (12), `count_upcrossings` counts samples with x[n] < x₀ ≤ x[n+1], and `variance_estimate_relative_sd` and `record_length_for_variance` are (13) and its inverse.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-19-random-processes-power-spectra`. The 21 tests:

1. **Analytic cases.** Both spectra integrate to σ² (`quad`, 10⁻⁸); Φ_w(0) = σ²L/(πV) and Φ_u(0) twice that; the ω⁻² tails have coefficients 3σ²V/(πL) and 2σ²V/(πL); R_w vanishes at 2L/V and has a negative lobe; the integral time scales are L/(2V) and L/V; the numerical Wiener–Khinchin transform of (6) returns (4)–(5) to 10⁻⁶; the filter gain squared equals πΦ to 10⁻¹⁰.
2. **Oracles.** `discretize_zoh` vs `cont2discrete`; `simulate_lti` vs `dlsim`; `welch_psd` vs `welch`; `periodogram` vs `periodogram`; `autocorrelation_estimate` vs `correlate`; the closed forms (11) vs `quad`: all to 10⁻¹⁰ or 10⁻⁸.
3. **Statistics of the seeded record** (each tolerance derived in a comment from (13) or its analogue): variance within 25 % (4.5σ), mean within 0.3 m/s (4σ), halves of the record agree, the fraction beyond 2σ is 4.5 % ± 4 %, Welch band averages within 30 %/30 %/10 % and the fitted high-frequency slope −2.0 ± 0.15, the periodogram's SD/mean in [0.7, 1.4] and Welch's below 0.4.
4. **Response.** σ_v within 15 % and σ_a within 8 % of (11) (the latter tight enough to catch the zero-order-hold bias); a is the derivative of v; a step gust gives a(0) = w/τ; counted zero and +1σ upcrossings within 25 % of Rice; the simulated exceedance fraction at 0.2 g within 25 % of erfc; the variance scatter over 16 seeds within a factor two of (13).
5. **Limits and errors.** τ → 0 and τ → ∞ in (11); monotonic dependence on τ; σ_a < σ_w/τ; `ValueError` on bad inputs.

### Compare

`main()` prints, for the lesson values, the band ratios Welch/Dryden of 0.989, 1.078 and 1.030 in the decades 0.01–0.1, 0.1–1 and 1–10 Hz, and the scatter SD/mean of the periodogram bins between 0.1 and 1 Hz, 1.03, against 0.22 for Welch with 23 segments of 8192 points (1/√23 = 0.21). `scipy.signal.welch` gives the same numbers to 14 digits, because it is the same algorithm; we wrote it to know what the window, the scaling and the segment count do, and we call the library in practice because it also handles odd lengths, other windows, two-sided and complex inputs, and `nfft` padding without our having to re-derive the bookkeeping each time. The same holds for `cont2discrete` and `dlsim` against `discretize_zoh` and `simulate_lti`: ours exist to show that "simulate a random process" means one matrix exponential and one loop, nothing more mysterious.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All records are synthetic, generated by `generate_gust` from seeded white noise; nothing here is a measurement.

1. **Scatter of the variance estimate.** For T = 100, 400 and 1600 s, generate 24 records (seeds 100–123, Δt = 0.1 s), compute each sample variance, and take the SD across seeds divided by σ_w². Before running, predict the three numbers from (13): 0.250, 0.125, 0.0625.
2. **Periodogram against Welch.** On the 2000 s record, compute both estimates, divide each by the Dryden model at every bin, and compute SD/mean of the ratio between 0.1 and 1 Hz. Predict: 1 for the periodogram, 1/√K for Welch.
3. **Resolution against scatter.** Repeat the Welch estimate with `nperseg` = 1024, 8192 and 65 536 and look at the lowest resolved frequency and the scatter. Which choice resolves the knee at 0.032 Hz, and what does it cost?
4. **Heave response.** Simulate v and a for the record, compute the Welch spectrum of Δn = a/g and overlay |H_a|²Φ_w/g²; count the zero crossings of v and compare with 0.148 s⁻¹ × 2000 s.

Raw data format: a CSV `T, predicted_rel_sd, measured_rel_sd` for item 1; the figure `gust_spectrum.png` (four panels: record, gust PSD with Welch and both Dryden forms, autocorrelation estimate against (6), load-factor PSD against |H_a|²Φ_w/g²) for the rest.

What you should find: the measured relative SDs are 0.236, 0.152 and 0.066 against the predicted 0.250, 0.125 and 0.0625, each within the ±15 % that an SD estimated from 24 samples carries (1/√(2 × 23)); the periodogram ratio has SD/mean 1.03 and Welch 0.22; and the load-factor spectrum peaks near √(bc)/(2π) = 0.12 Hz, where the rising |H_a|² meets the falling Φ_w: the aircraft turns a gust spectrum that is all below 0.03 Hz into bumps at 0.1–0.3 Hz.

## 5. Validate and improve

| quantity (lesson aircraft, seed 2019, 2000 s) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| ∫Φ_w dω, ∫Φ_u dω | σ² = 2.25 | 2.250000 | 0 | exact forms |
| Φ_w(0) | 3.581 | 3.580986 | 0 | exact |
| sample SD of the record (m/s) | 1.5 ± 0.042 (2.8 %) | 1.5077 | +0.5 % | 0.2σ scatter |
| sample mean (m/s) | 0 ± 0.075 | −0.060 | 0.8σ | scatter |
| Welch/Dryden, three decades | 1 ± 0.07, 0.05, 0.01 | 0.989, 1.078, 1.030 | ≤ 1.5σ | segment scatter |
| periodogram SD/mean | 1 | 1.03 | 3 % | 1800 bins |
| τ (s) | 0.3501 | 0.350140 | 0 | exact |
| σ_v (m/s) | 1.426 | 1.4262 closed form, 1.4302 simulated | +0.3 % | record scatter |
| σ_a (m/s²), σ_n (g) | 1.328, 0.135 | 1.3275 closed form, 1.3630 simulated | +2.7 % | 2σ of the a-variance scatter (1.3 % in SD); seeds 0–5 average −0.4 % |
| fraction beyond 0.3 g | 2.67 % | 3.09 % simulated | +16 % | the +2.7 % in σ_a moves a 2.2σ tail by 16 % |
| zero upcrossings of v in 2000 s | 296 (Rice) | 301 | +1.7 % | √296 = 17 scatter |
| record for 2 % variance (min) | 260 | 260.4 | 0 | (13) |

Every difference has a named cause and a size consistent with the scatter derived in 2.8. Two of them deserve a second look. The simulated σ_a is 2.7 % high on this seed; over seeds 0–5 the ratio is 0.991 ± 0.025, so the method is unbiased and this seed is a two-sigma draw. Had `simulate_heave` used a zero-order hold, the ratio would have been 1.056 on every seed: a *systematic* 5 %, visible only because the closed form (11) exists. That is the pattern to carry forward: a Monte Carlo result without an independent closed form or limiting case cannot tell a bias from a fluctuation. The exceedance fraction at 0.3 g inherits the σ_a error amplified by the tail: a 2.2σ level is where a few per cent in σ become 16 % in probability, and in a certification analysis that sensitivity is the reason the gust intensity is treated as uncertain rather than fixed.

What the agreement does not show: that real turbulence is Dryden. It shows that the generator produces what the Dryden model specifies and that the estimators recover it. Against flight data the ω⁻² tail would be the first thing to fail (von Kármán's −5/3 fits measured spectra better), the Gaussian tails the second, and stationarity the third. The record length (13) is the quantity to take from this lesson to a real test: 260 minutes for 2 % is not available, so a flight-test spectrum is a 10–20 % measurement and should be reported as one.

Improvements worth making: implement the von Kármán spectrum and generate it through a rational approximation (extension 1); add pitch dynamics so that σ_n includes the short-period response; replace the frozen scale length by the specification's altitude dependence. Improvement worth *not* making: refining Δt below 0.02 s, which changes σ_a by 0.1 % while the statistical scatter is 1 %.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of the normalisation of (5) and of σ_v² in (11) with its two limits; `gust_spectrum.png`; the scatter table; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the RMS load factor, the exceedance fraction at 0.3 g, the zero-crossing rate of v, the record length for 2 % and the assumption you trust least.

### Adversarial self-checks

- **What if σ_w doubles?** Every spectrum and variance quadruples, every RMS doubles, the crossing rate ν₀ is unchanged (it is a ratio of RMS values), and the fraction beyond 0.3 g rises from 2.7 % to 27 %. What if L doubles? Φ(0) doubles and b halves; σ_v rises slightly (the aircraft follows more of the gust) and σ_a falls by ≈ √2 when r ≫ 1. Check both with the code.
- **What if V → 0? V → ∞?** b → 0 and c → 0 together; σ_a² ≈ (3/2)σ_w²bc → 0 as V². At high V the quasi-steady and plunge-only assumptions fail long before the formula does.
- **Which assumption fails first in the real world?** Stationarity of the 33-minute record; then the ω⁻² tail if anything above 1 Hz (a structural mode, an accelerometer's crossing rate) is asked of the model.
- **What measurement could falsify the model?** A gust-probe record with a well-resolved spectrum between 0.1 and 10 Hz: a slope of −5/3 instead of −2 falsifies Dryden in favour of von Kármán, and a non-zero skewness or an excess of 4σ events falsifies the Gaussian assumption.
- **Does the result survive a change in Δt, record length or seed?** σ_v, σ_a and the band ratios move within the derived scatter for seeds 0–5 and for Δt = 0.005 s; the zero-order-hold variant of `simulate_heave` does not survive the comparison with (11). Say which numbers in your note you would quote to two figures and which to one.
- **Which decision would change if the uncertainty doubled?** If σ_w were known only to ±40 % instead of ±20 % (0.6 % to 6.5 %), the exceedance fraction at 0.3 g would range from 0.02 % to 11 %: the "acceptable ride" decision at 0.1 g RMS (which needs V ≤ 44 m/s, from σ_n ∝ V) would be replaced by "slow down and accept a wide band", and the record-length formula says that halving that uncertainty needs four times the data.

## Check yourself

1. Define wide-sense stationarity in one sentence and name the two properties of R(τ) that follow from the definition alone.
2. State the Wiener–Khinchin relation in this lesson's one-sided convention and the check that catches a wrong factor of 2 or π.
3. Write the Dryden vertical-gust spectrum and give its value at ω = 0 and its behaviour as ω → ∞ in terms of σ_w, L and V.
4. Why does sampled white noise need variance 1/Δt, and what happens to the generated gust's variance if you use variance 1 instead?
5. The periodogram of a 2000 s record and of a 20 000 s record have the same scatter. Why, and what does Welch's method trade to reduce it?
6. For the lesson aircraft, why is the RMS load factor (0.135 g) so much smaller than the sharp-edged-gust value (0.437 g)?
7. Rice's formula gives a finite zero-crossing rate for v but not for w_g or a. What property of the Dryden spectrum is responsible?
8. How long a record does a 5 % relative error in σ̂_w² need for L = 300 m and V = 60 m/s, and how does the answer scale with V?
