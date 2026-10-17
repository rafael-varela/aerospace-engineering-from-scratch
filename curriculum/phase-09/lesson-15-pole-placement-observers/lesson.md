# P09.15  Pole Placement & Observers

> If you can measure every state, you can put every closed-loop pole where you like. If you cannot, build a model that estimates the states and feed back the estimates. The two designs do not interfere.

| | |
|---|---|
| **Learning outcomes** | 1. Derive Ackermann's formula K = e_nᵀW_c⁻¹φ_d(A) for single-input full-state feedback, the reference gain N for zero steady-state error, and the Luenberger observer gain by duality, L = ackermann(Aᵀ, Cᵀ)ᵀ. 2. Prove the separation principle: the loop with u = Nr − Kx̂ has exactly the regulator poles eig(A − BK) and the observer poles eig(A − LC). 3. Implement and test pole placement and observers from scratch, verify eig(A − BK) and eig(A − LC) against their targets and against `scipy.signal.place_poles`, and show by fitting that the observer error decays at the slowest placed observer pole. 4. Choose regulator and observer speeds for a reaction-wheel satellite with a failed gyro by trading wheel torque against settling, and estimation speed against SYNTHETIC star-tracker noise. |
| **Prerequisites** | `P09.14` (controllability and observability, PBH, duality), `P09.13` (state-space models, similarity transforms, Cayley–Hamilton), `P09.04` (ζ, ω_n and overshoot), `P02.05` (matrix exponential and zero-order hold), `P02.19` (random processes: white noise, variance). Skill check: you can write the characteristic polynomial of a 2 × 2 matrix and match its coefficients to s² + 2ζω_n s + ω_n². |
| **Engineering question** | A small satellite's gyro has failed; the star tracker still measures the attitude angle about one axis. Can a reaction-wheel controller still slew 10° and hold, with the rate supplied by an observer, without exceeding a 10 mN m wheel limit? How fast should the observer be? Answer, for a representative 3 kg m² axis: regulator poles at ζ = 0.707, ω_n = 0.1 rad/s need K = (0.030, 0.424) and a peak torque of 5.24 mN m; observer poles at −0.3 and −0.6 s⁻¹ give L = (0.9, 0.18), and the estimation error decays at a fitted 0.29998 s⁻¹ against 0.3 placed. Doubling the observer speed roughly triples the rate-estimate noise: 0.167 → 0.485 → 1.29 → 3.47 arcsec/s for slowest poles 0.15, 0.3, 0.6, 1.2 s⁻¹ with 10 arcsec SYNTHETIC tracker noise. |
| **Assumptions and limits** | Linear, time-invariant, single-input plants with full knowledge of A, B, C. The satellite axis is a rigid double integrator Jθ̈ = τ (no flexible modes, no gyroscopic coupling, no disturbance torques), the star tracker is an angle measurement with white SYNTHETIC noise held at 10 Hz, and the wheel is linear below its limit. The four-state pitch example is the representative model of `P09.11`. Ackermann's formula inverts the controllability matrix, which becomes ill-conditioned for high-order or weakly controllable systems; robust MIMO placement (as in `place_poles`) and LQR (`P09.16`) are the remedies. |
| **Analytical work** | The closed-loop characteristic polynomial of A − BK for the double integrator, by hand; Ackermann's formula through the controllable canonical form and Cayley–Hamilton; the uniqueness of K for a single input; duality for L; the reference gain N; the separation principle from the (x, e) coordinates; the error dynamics ė = (A − LC)e − Lv and its stationary covariance; the peak torque Nr = Jω_n²r at t = 0. |
| **Computational deliverable** | `starter/main.py`: `char_poly`, `ackermann`, `observer_gain`, `reference_gain`, `observer_closed_loop`, `simulate_observer_loop`, `fit_decay_rate` (`expm`, `discretize_zoh`, `satellite_axis`, `pitch_model`, `second_order_poles`, `main` are given). Tests prove: **eig(A − BK) equals the target poles** for five systems (satellite, four-state pitch, random 3- and 5-state); **K equals `scipy.signal.place_poles`** (unique for a single input); the satellite K is J(ω_n², 2ζω_n); placing the open-loop poles gives K = 0; K transforms as KT under x = Tz; uncontrollable and unobservable cases raise; eig(A − LC) equals its targets and L equals `place_poles` on the dual; the 8-state loop's eigenvalues are the union of the two sets (separation); N removes the steady error; the observer error equals e^{(A−LC)t}e₀ whatever K and r are; a perfect initial estimate reproduces state feedback; **the observer error's fitted decay rate equals the slowest placed observer pole to 10⁻⁴** for three pole sets; peak torque grows as ω_n²; a faster observer passes more SYNTHETIC tracker noise. |
| **Practical activity** | **computational experiment**: regulator-speed versus wheel-torque sweep, observer-speed versus SYNTHETIC tracker-noise sweep (with a covariance prediction to compare), and a fitted-decay check of the observer. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `pole_placement.png`, the derivations of Ackermann's formula, N and the separation principle, the two sweep tables (the noise one labelled SYNTHETIC), and a half-page recommendation of regulator and observer poles for the failed-gyro satellite. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: Ackermann's formula and separation are derived; all public tests pass; the recommendation states the peak torque, the settling behaviour, the fitted observer decay rate and the SYNTHETIC noise level of the chosen design, and why the observer is not made faster. |
| **References** | `@mit-ocw-16-30-feedback-control-systems` Topic 11 "Full-state feedback control" (https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec11/): Full-state Feedback Controller; Ackermann's Formula; Origins of Ackermann's Formula; Reference Inputs; Pole Placement Examples (required-free): coefficient matching, K = [0 … 1]W_c⁻¹φ_d(A) and the reference gain. Topic 14 "Open-loop and closed-loop estimators" (.../mit16_30f10_lec14/): Closed-loop Estimator; Estimation Gain Selection; Dual Design Approach; Where Put Estimator Poles? (required-free): the observer gain by duality, ė = (A − LC)e, and observer poles 2–5 times faster than the regulator's. Topic 15 "Combined estimators and regulators" (.../mit16_30f10_lec15/): Combined Estimators and Regulators; The Compensator; Satellite Control Design (required-free): the separation of the closed-loop poles, on a satellite example. `@mit-ocw-6-011-signals-systems-and-inference` Lecture 9 "Observers for State Estimation"; Lecture 10 "Observers, State Feedback" (https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/resources/mit6_011s18lec10/); Lecture 11 "State Feedback, Observer-Based Feedback" (supplementary-free; sparse slides): the same results in a second notation. `@mit-ocw-18-06-linear-algebra` Lecture 21 "Eigenvalues and eigenvectors", Lecture 23 "Differential equations and exp(At)" (supplementary-free review). `@scipy-userguide` API reference "Signal processing (scipy.signal)": `place_poles` (supplementary-free): the test oracle only. |
| **Extensions** | ★ Add a constant disturbance torque (solar pressure) and an integral state to the regulator; place the extra pole and show zero steady error. ★ Build a reduced-order (one-state) observer that uses the measured angle directly and estimates only the rate; compare its noise with the full observer. ★ Place the poles of the 4-state pitch loop with `place_poles` for a two-input version (elevator and throttle) and compare its gain norm with the single-input Ackermann gain. |
| **Time** | unmeasured |

## 1. The engineering problem

`P09.11` and `P09.12` designed controllers that see only the output and shape a single loop. When the whole state is available, there is a more direct method: feed back every state, u = −Kx, and the closed-loop matrix becomes A − BK. If (A, B) is controllable (`P09.14`), its eigenvalues can be put anywhere you choose. The gain needs only a matrix formula. That is **pole placement**.

The catch is "the whole state". A satellite's attitude controller wants angle and rate, but after a gyro failure only the star tracker's angle remains. Differentiating the noisy angle to get a rate amplifies the noise without limit. The alternative is an **observer**: a copy of the plant model, driven by the same commands, whose state is nudged toward the measurement by a gain L. Its error obeys its own linear equation, with eigenvalues eig(A − LC) that we place exactly like the regulator's, by **duality**. The **separation principle** then says the two designs do not interfere: the closed loop has the regulator poles and the observer poles, nothing else. What does interfere is noise. A fast observer trusts the measurement more and passes more of its noise into the rate estimate, and from there into the wheel.

Success criteria, written first: eig(A − BK) and eig(A − LC) equal their targets to better than 10⁻⁶ relative and agree with `scipy.signal.place_poles`; the 4-state loop with the observer has exactly the four placed poles; a 10° slew stays within the 10 mN m wheel limit and has a few percent overshoot; the observer error, fitted on a log scale, decays at the slowest placed observer pole to 10⁻⁴.

## 2. Understand and derive

### 2.1 Full-state feedback on a double integrator, by hand

One satellite axis, x = (θ, ω), Jθ̈ = τ: A = [[0, 1], [0, 0]], B = (0, 1/J)ᵀ. With τ = −k₁θ − k₂ω, A − BK = [[0, 1], [−k₁/J, −k₂/J]], whose characteristic polynomial is s² + (k₂/J)s + k₁/J. Matching s² + 2ζω_n s + ω_n²:

**k₁ = Jω_n², k₂ = 2ζω_nJ.**   (1)

This is a PD controller with gains chosen by pole location. For J = 3 kg m², ζ = 0.707 and ω_n = 0.1 rad/s: **K = (0.030 N m/rad, 0.4242 N m s/rad)**.

### 2.2 Ackermann's formula

For a general single-input (A, B), transform to the controllable canonical form of `P09.13`. In those coordinates A − BK has the last row of A replaced by (open-loop coefficients minus K), so each gain entry moves one coefficient of the characteristic polynomial: placing poles is coefficient matching. Transforming back gives a coordinate-free result. With W_c = [B, AB, …, A^{n−1}B] and the desired polynomial φ_d(s) = sⁿ + a₁sⁿ⁻¹ + ⋯ + a_n,

**K = e_nᵀ W_c⁻¹ φ_d(A),**   (2)

(16.30 Topic 11, "Ackermann's Formula" and "Origins of Ackermann's Formula"), where φ_d(A) = Aⁿ + a₁Aⁿ⁻¹ + ⋯ + a_nI and e_nᵀ = (0, …, 0, 1). A direct check: by Cayley–Hamilton applied to A − BK, φ_d(A − BK) = 0. Expanding φ_d(A − BK) and multiplying by e_nᵀW_c⁻¹ leaves φ_d(A) − K-terms that vanish only for (2). The formula needs W_c⁻¹: **pole placement is possible exactly when (A, B) is controllable**, which is the `P09.14` test, and `ackermann` raises otherwise. For a single input the n entries of K are fixed by the n coefficients, so **K is unique**. This is why the test can demand equality with `scipy.signal.place_poles`, which uses a completely different algorithm (orthogonal transformations, chosen for robustness in the multi-input case).

The **reference gain**: with u = Nr − Kx, the steady state of ẋ = (A − BK)x + BNr is x_ss = −(A − BK)⁻¹BNr, and y_ss = r requires

**N = −1/(C(A − BK)⁻¹B).**   (3)

For the satellite, N = k₁ = 0.030: the reference simply replaces θ in the proportional term.

### 2.3 The Luenberger observer by duality

The observer is a copy of the model corrected by the output error:

x̂̇ = Ax̂ + Bu + L(y − Cx̂).   (4)

Subtracting from ẋ = Ax + Bu, the error e = x − x̂ obeys **ė = (A − LC)e**, independent of u. Its eigenvalues are those of (A − LC)ᵀ = Aᵀ − CᵀLᵀ: a pole-placement problem for the dual pair (Aᵀ, Cᵀ) (16.30 Topic 14, "Dual Design Approach"). So **L = ackermann(Aᵀ, Cᵀ, observer poles)ᵀ**, which exists exactly when (A, C) is observable. For the star-tracker-only satellite (C = [1, 0]) with observer poles −p₁, −p₂, coefficient matching gives L = (p₁ + p₂, p₁p₂): for −0.3 and −0.6, **L = (0.9, 0.18)**. With a gyro only (C = [0, 1]) the pair is unobservable and no L exists (`P09.14`), which is why losing the star tracker would be worse than losing the gyro.

### 2.4 The separation principle

Close the loop on the estimate, u = Nr − Kx̂, and use (x, e) as coordinates (x̂ = x − e):

```
d/dt [x]   [A − BK    BK  ] [x]   [BN]
     [e] = [  0     A − LC] [e] + [ 0] r.    (5)
```

The matrix is block upper triangular, so its eigenvalues are **eig(A − BK) ∪ eig(A − LC)**. The regulator and observer can be designed separately (16.30 Topic 15, "Combined Estimators and Regulators"). Three consequences:

- the observer error never depends on r or K (it is driven by nothing but its initial value and noise);
- the transfer function from r to y is the same as with full-state feedback (the observer modes are uncontrollable from r: they cancel);
- if x̂(0) = x(0), the response is exactly the state-feedback response.

For the satellite with the gains above, the 4 × 4 loop has eigenvalues −0.0707 ± 0.0707j, −0.3 and −0.6, as `main()` prints.

### 2.5 Noise: why the observer is not made arbitrarily fast

With measurement noise v, y = Cx + v and (4) gives **ė = (A − LC)e − Lv**. Larger L means faster decay of an initial error but more noise injected. For white noise of spectral density S, the stationary error covariance P solves the Lyapunov equation (`P09.14`)

(A − LC)P + P(A − LC)ᵀ + LLᵀS = 0.   (6)

For the double integrator with poles −ω₀, −2ω₀, L = (3ω₀, 2ω₀²), and (6) gives a rate-error standard deviation proportional to ω₀^{3/2}. Each doubling of the observer speed multiplies the rate noise by 2^{1.5} = 2.8. A star tracker with σ = 10 arcsec sampled and held at 10 Hz has S ≈ σ²Δt. The predicted rate-error standard deviations are 0.150, 0.424, 1.200 and 3.394 arcsec/s for ω₀ = 0.15, 0.3, 0.6 and 1.2 s⁻¹ (a probe solving (6) with `scipy.linalg.solve_continuous_lyapunov`; §4 compares them with the SYNTHETIC simulation).

### 2.6 Where the model stops being valid

- **Model error.** The observer uses the model's A and B. If the true inertia is 20 % larger, the loop poles move: with true J = 3.6 kg m² the eigenvalues of the actual loop are −0.552, −0.386 and −0.052 ± 0.066j instead of the placed set (Problem 5). Separation is a property of the model, not of the hardware.
- **Disturbances.** A constant disturbance torque (solar pressure, a residual magnetic dipole) biases both the attitude and the estimate. Pure state feedback with N has no integral action; extension 1 adds it.
- **Torque limit.** The peak torque of a step is Nr at t = 0. It grows as ω_n², so fast regulator poles saturate the wheel, and then the linear design no longer describes the response.
- **Flexible modes.** Solar arrays (`P09.13` Problem 4) add lightly damped modes that the rigid model does not have; an observer tuned too fast can excite them.

## 3. Build the computational model

Open `starter/main.py`. `expm`, `discretize_zoh`, `satellite_axis`, `pitch_model` and `second_order_poles` are given.

### Step 1: placement

```python
def ackermann(A, B, poles):
    n = A.shape[0]
    Wc = np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(n)])
    # raise ValueError if the smallest singular value of Wc is ~0
    phi = np.zeros((n, n))
    for c in char_poly(poles):            # Horner evaluation of phi_d(A)
        phi = phi @ A + c * np.eye(n)
    row = np.linalg.solve(Wc.T, np.eye(n)[-1])   # e_n^T Wc^{-1} without an explicit inverse
    return (row @ phi)[None, :]
```

`char_poly` is `np.poly` with a check that the coefficients are real (complex poles must come in pairs). `observer_gain` is `ackermann(A.T, C.T, poles).T`, re-raising with an observability message. `reference_gain` is (3).

### Step 2: the loop with the observer

`observer_closed_loop` builds the 2n × 2n matrix in (x, x̂) coordinates, [[A, −BK], [LC, A − BK − LC]], with inputs (r, v): (BN; BN) for the reference and (0; L) for the noise. `simulate_observer_loop` discretises it once by ZOH and steps it; with r and v held over each sample this is exact. It returns x, x̂, u = Nr − Kx̂ and e = x − x̂.

### Step 3: the fitted decay

`fit_decay_rate(t, e, t_start, t_end)` fits log‖e‖ = c − σt by least squares over the window. The window must start after the faster observer mode has died out (here (fast − slow) × t ≥ 12 or so) and end before ‖e‖ reaches the rounding floor of the states it is the difference of. A window that runs into that floor bends the fit: an early version of `main()` fitted from 30 to 120 s and got 0.272 instead of 0.300.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-15-pole-placement-observers`. The tests:

1. **Independent check, part 1**: eig(A − BK) equals the targets for five systems (two satellite designs, the 4-state pitch model with poles −4 ± 4j, −6, −20, random 3- and 5-state systems); K equals `scipy.signal.place_poles(...).gain_matrix` for four of them to 10⁻⁷.
2. **Structure**: K = J(ω_n², 2ζω_n) for the satellite; targets equal to the open-loop eigenvalues give K = 0; under x = Tz the gain becomes KT; an uncontrollable pair, a gyro-only observer and an incomplete complex pair all raise.
3. **Observer**: eig(A − LC) equals its targets and L equals `place_poles(Aᵀ, Cᵀ)` transposed, for the satellite and the pitch model.
4. **Separation**: the 8 × 8 pitch loop has exactly the union of the regulator and observer poles; N gives θ → r; the error trajectory is identical for two different (K, r) pairs and equals e^{(A−LC)t}e₀; a perfect initial estimate gives zero error and the state-feedback trajectory.
5. **Independent check, part 2**: for observer poles (−0.3, −0.6), (−0.5, −1.5) and (−0.2, −0.45) the fitted decay rate equals the slowest pole to 10⁻⁴; the fit is exact on known exponentials.
6. **Trade-offs**: peak torque Jω_n²r for ω_n = 0.05, 0.1, 0.2; with seeded SYNTHETIC noise, an observer four times faster has more than three times the rate-estimate noise.

### Compare

`scipy.signal.place_poles` is the oracle. For a single input it must return our K exactly, and it does to 10⁻⁷ even for the 4-state pitch model, where W_c has entries of very different sizes. For several inputs K is not unique, and `place_poles` picks the one that makes the closed-loop eigenvectors as well conditioned as possible, which Ackermann cannot do. `main()` prints:

```
regulator poles [-0.0707+0.0707j -0.0707-0.0707j] -> K = [0.03   0.4242] (= J (wn^2, 2 zeta wn))
observer poles [-0.3 -0.6] -> L = [0.9  0.18]; reference gain N = 0.03000 (= K[0])
eig of the 4x4 loop with the observer: [-0.6   +0.j       -0.3   +0.j       -0.0707-0.070721j -0.0707+0.070721j]
10 deg slew from a 0.05 deg/s tip-off: overshoot 4.46 %, peak torque 5.236 mN m, final error 2.08e-04 deg
same slew with a perfect initial estimate: peak torque 4.866 mN m; max |theta_observer - theta_state_feedback| = 0.1618 deg
observer error: fitted decay rate 0.29998 1/s over 20-60 s (slowest placed observer pole 0.30); |e| falls from 8.73e-04 to 5.17e-11 at 60 s
```

The observer starts with x̂ = 0 while the satellite has an unknown 0.05°/s tip-off rate (a residual rate after separation from the launcher, a representative value). The estimate converges within about 15 s. The attitude differs from the perfect-knowledge response by at most 0.16°, and the overshoot is 4.46 %, close to the 4.3 % of ζ = 0.707.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment**, with a SYNTHETIC noise record. Write each prediction down first.

1. **Regulator speed against wheel torque.** For ω_n = 0.05, 0.1, 0.2, 0.4 rad/s (ζ = 0.707), simulate a 10° step with a perfect estimate and record the peak torque. Predict it from Nr = Jω_n²r and mark which designs exceed 10 mN m.
2. **Observer speed against noise.** With the regulator fixed (ω_n = 0.1 rad/s) and observer poles (−ω₀, −2ω₀) for ω₀ = 0.15, 0.3, 0.6, 1.2 s⁻¹, hold the commanded attitude at 0 for 120 s with SYNTHETIC tracker noise (σ = 10 arcsec, white, held at 10 Hz, seed 1509). Record the standard deviation of the rate-estimate error and of the wheel torque after 40 s. Predict the rate error from (6) and the ω₀^{3/2} law.
3. **Decay check.** Start the observer at x̂ = 0 with a 0.05°/s tip-off and fit the decay of ‖e‖ over 20–60 s.

Raw data format: CSV tables `wn_rad_s, peak_torque_mNm, predicted_mNm`; `omega0_per_s, rate_err_arcsec_s, rate_err_predicted, torque_rms_uNm` (header: SYNTHETIC star-tracker noise, σ = 10 arcsec, 10 Hz, seed 1509); `window_s, fitted_rate, placed_rate`.

What you should find: peak torques 1.309, 5.236, 20.94 and 83.78 mN m, equal to Jω_n²r: only ω_n ≤ 0.1 rad/s stays within the wheel limit (a faster design would saturate, and the linear analysis would no longer hold). SYNTHETIC rate-estimate errors 0.167, 0.485, 1.286 and 3.474 arcsec/s against 0.150, 0.424, 1.200 and 3.394 predicted by (6). The ratio per doubling is 2.9, 2.65 and 2.7 (2.8 predicted). The simulation sits 2–14 % above the continuous prediction, because the noise is held over each 0.1 s sample, which is only approximately white at these observer speeds, and the statistics come from 800 correlated samples. The torque noise grows from 0.586 to 7.78 μN m. The fitted decay rate is 0.29998 s⁻¹. An observer at (−0.3, −0.6), three to six times faster than the regulator, settles its error long before the slew is done and adds under 1.4 μN m of torque noise. Faster buys nothing the regulator can use.

## 5. Validate and improve

| quantity | hand (§2) | code | independent check | difference |
|---|---|---|---|---|
| satellite K | J(ω_n², 2ζω_n) = (0.030, 0.4242) | (0.03, 0.4242) | `place_poles` (test) | < 10⁻⁷ |
| satellite L | (p₁ + p₂, p₁p₂) = (0.9, 0.18) | (0.9, 0.18) | `place_poles` on the dual (test) | — |
| loop eigenvalues | union, §2.4 | −0.0707 ± 0.0707j, −0.3, −0.6 | 8-state pitch union (test) | < 10⁻⁵ |
| observer decay rate | 0.3 s⁻¹ | fitted 0.29998 s⁻¹ | three pole sets to 10⁻⁴ (test) | fast-mode residue |
| peak torque, 10° step | Jω_n²r = 5.236 mN m | 5.236 mN m | — | — |
| rate-error noise, ω₀ = 0.3 | (6): 0.424 arcsec/s | SYNTHETIC 0.485 arcsec/s | — | held noise, finite record |

The placement and the separation principle are confirmed to rounding by three routes: hand coefficient matching, Ackermann's formula and SciPy's independent algorithm. The decay rate is confirmed by fitting the simulated error, a measurement of the design rather than a restatement of it. The noise result is a SYNTHETIC simulation compared with a covariance prediction, which agree to within 15 %, not a measurement of a real star tracker. What it shows is the trade-off, not the exact number for any particular sensor.

Improvements worth making: integral action for disturbance torques (extension 1); a reduced-order observer (extension 2); choosing the observer gain from the noise statistics rather than from pole locations, which is the Kalman filter and beyond this phase. Improvement worth not making: pushing the regulator poles faster to "tighten" the slew; the wheel saturates beyond ω_n ≈ 0.14 rad/s for 10° steps.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2), (3), L by duality and the separation principle (5); `pole_placement.png`; the three experiment tables, the noise table labelled SYNTHETIC; and a half-page recommendation (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) of regulator and observer poles for the failed-gyro satellite: K, L, N, the peak torque, the overshoot, the fitted observer decay and the SYNTHETIC rate noise, with one sentence on what you would check first in flight.

### Adversarial self-checks

- **What if the inertia is 20 % larger than the model?** The actual loop's poles move to −0.552, −0.386 and −0.052 ± 0.066j (Problem 5): slower and less damped regulator modes, but still stable, and θ still settles on r because the plant integrates. Separation holds only for the model.
- **What if the step is 30° instead of 10°?** The peak torque Nr = 15.7 mN m exceeds the limit. Either slow the poles, shape the reference (a ramp instead of a step), or accept saturation and add anti-windup thinking from `P09.11`.
- **Which assumption fails first in the real world?** A disturbance torque: without integral action it leaves a constant attitude error of τ_d/k₁, which for k₁ = 0.03 N m/rad means 1.9° per mN m.
- **What measurement could falsify the design?** Log the observer's innovation y − Cx̂ in flight. If its mean is not zero or its spectrum is not flat, the model (J, a bias, a flexible mode) is wrong.
- **Does the result survive a change of algorithm?** Ackermann and `place_poles` agree to 10⁻⁷; the simulation by ZOH is exact for held inputs; the fitted decay matches the placed pole.
- **Which decision changes if the tracker noise doubles?** The rate noise doubles at every observer speed; to keep the same rate noise the observer must be slowed by 2^{2/3} = 1.6, which still leaves it twice as fast as the regulator for the (−0.3, −0.6) choice. If the noise rose tenfold, the observer would have to approach the regulator's speed and separation would stop being a comfortable approximation in practice.

## Check yourself

1. Write Ackermann's formula and say what property of (A, B) it requires.
2. For Jθ̈ = τ with full-state feedback, give k₁ and k₂ for poles with damping ζ and natural frequency ω_n.
3. Why is the single-input gain K unique, and why is the multi-input gain not?
4. How is the observer gain L obtained from a pole-placement routine, and what property of (A, C) does it require?
5. State the separation principle and the coordinate change that proves it.
6. Why does the observer error not depend on the reference or on the regulator gain?
7. An observer's poles are −0.3 and −0.6 s⁻¹. At what rate does the estimation error decay after the first few seconds, and how would you measure it?
8. Why not make the observer ten times faster than the regulator?
