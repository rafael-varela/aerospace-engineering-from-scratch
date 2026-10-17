# P09.18  Robustness & Uncertainty ★

> Every model is wrong at high frequency. Robustness analysis asks how wrong it may be before the loop goes unstable, answers with one number, and tells you honestly when that number is too cautious.

| | |
|---|---|
| **Learning outcomes** | 1. Define the sensitivity functions S, T, KS and GS of a feedback loop, state what each one governs (disturbance rejection, noise, control effort, robustness), and derive S + T = 1. 2. Compute the H∞ norm (the peak gain over frequency) of a stable system two independent ways, by a refined frequency sweep and by bisection on the imaginary-axis eigenvalues of a Hamiltonian matrix, and use it to obtain Ms and Mt. 3. Derive the disk-margin guarantees from Ms (gain in (Ms/(Ms + 1), Ms/(Ms − 1)), phase ≥ 2 arcsin(1/(2Ms))) and check them against measured margins. 4. Model unmodelled actuator dynamics as multiplicative uncertainty G = G₀(1 + WΔ), apply the small-gain test ‖WT‖∞ < 1, construct the destabilising perturbation when it fails, and judge how conservative the test is against the true plant. 5. Verify the Bode sensitivity integral ∫₀^∞ ln|S| dω = πΣp_i numerically and use it to state what no controller can achieve. |
| **Prerequisites** | `P09.10` (Nyquist criterion, gain, phase and delay margins), `P09.12` (lead compensation of a servo), `P09.08` (sensitivity and complementary sensitivity), `P09.13` (state space, transfer function to state space), `P09.05` (frequency response), `P01.16` (eigenvalues), `P02.07` (complex numbers). Skill check: you can sketch the Nyquist plot of K/(s(τs + 1)) and read off its phase margin. |
| **Engineering question** | A control-surface servo is designed on a two-parameter motor model, but the real actuator has a gear-compliance mode near 150 rad/s and a 5 ms transport delay that the model omits. For a family of lead compensators of increasing gain, which ones are guaranteed to stay stable on any actuator within a stated uncertainty bound, which are stable on the real actuator, what margins does each design's sensitivity peak guarantee, and how large is the gap between the guarantee and the truth? |
| **Assumptions and limits** | Linear time-invariant single-input single-output loop with negative unity feedback; the nominal servo G₀ = 30/(s(0.04s + 1)) (a representative small actuator, shaft angle per volt); the "true" actuator is G₀ times a second-order mode (150 rad/s, ζ = 0.3) times e^{−0.005s}, approximated by a Padé rational function where a polynomial is needed; the uncertainty is unstructured (any stable complex Δ with ‖Δ‖∞ ≤ 1), which makes the small-gain test conservative for real-parameter errors; no saturation or other nonlinearity (`P09.19`). The Bode integral needs a loop of relative degree at least two. |
| **Analytical work** | S, T, KS, GS from the loop equations; S + T = 1; the H∞ norm as a peak and its Hamiltonian characterisation; the disk of radius 1/Ms around −1 and the margins it implies; multiplicative uncertainty, the factorisation 1 + CG = (1 + L₀)(1 + WTΔ) and the small-gain theorem, with the all-pass construction showing necessity; the Bode sensitivity integral and the waterbed bound ln Ms ≥ (πp − ω₁ ln ε)/(ω₂ − ω₁). |
| **Computational deliverable** | `starter/main.py`: `closed_loop_tfs`, `is_closed_loop_stable`, `tf_to_ss`, `hinf_norm_sweep`, `hamiltonian_has_imaginary_eigs`, `hinf_norm_bisection`, `disk_margin_bounds`, `margins`, `robust_stability_index`, `destabilising_perturbation`, `bode_sensitivity_integral` (`freqresp`, `pade`, the plant and weight models, `lead`, `tf_hinf_norm`, `weight_covers`, `unstable_airframe`, `design_table` and `main` are given). Tests prove: the canonical realisation reproduces `scipy.signal.freqs`; S + T = 1, S = 1/(1 + L), KS and GS pointwise, and the DC and high-frequency limits; the closed-loop stability test on four loops; the H∞ norm by sweep equals the Hamiltonian bisection to 10⁻⁶ on random systems and the closed-form resonance peak 1/(2ζ√(1 − ζ²)) to 10⁻⁷; the Hamiltonian test switches at the norm; measured margins of a textbook loop (GM = 3 exactly); the disk-margin formula and its Ms → 1 limit; measured margins satisfy the Ms bounds on five loops, including the true servo and an unstable airframe with a gain-reduction margin of 1/6; the weight covers the true actuator error; ‖WT‖∞ < 1 implies stability for 20 random covered perturbations; ‖WT‖∞ > 1 yields an all-pass Δ of size below 1 that puts a pole on the imaginary axis; the Bode integral is 0 for the stable servo and 2π for the airframe with two different controllers. |
| **Practical activity** | **computational experiment**: a design sweep of six lead compensators recording Ms, Mt, ‖WT‖∞, the disk bounds and the margins of the true actuator; the conservatism of the small-gain test (robust gain limit against true gain limit, Problem 6); the waterbed of an unstable airframe. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `robustness.png`, the derivations of S + T = 1, the disk margin and the small-gain factorisation, the design table, and a half-page note choosing a lead gain with its guaranteed and true margins and the evidence you would want from a hardware test. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the disk margin and small-gain results are derived; all public tests pass; the design table is explained, including a design that fails the test but is stable and one that passes; the note separates guarantees from measurements. |
| **References** | See part 2 and the list at the end of the lesson. |
| **Extensions** | ★ Replace the first-order weight by a second-order one shaped around the compliance mode and recompute the robust gain limit; how much of the conservatism of Problem 6 disappears? ★ Add a real gain uncertainty k ∈ [0.8, 1.2] and compute the exact real stability margin by scanning k, against the complex small-gain bound. ★ Compute the structured singular value (μ) lower bound for two simultaneous perturbations (gain and delay) by a power iteration. |
| **Time** | unmeasured |

## 1. The engineering problem

A servo that positions an elevator, a fin or a thrust-vector gimbal is usually designed on a model with two parameters: a motor gain and a mechanical time constant, G₀(s) = K_m/(s(τ_ms + 1)), identified from a step test. The lead compensator of `P09.12` gives that model a comfortable phase margin. Then the servo is built, and the real actuator has more: a gear train and linkage that act as a spring (a lightly damped mode), amplifier dynamics, a sampling and computation delay. None of these appears in G₀, and all of them matter near and above the crossover frequency. A design pushed to a high crossover for speed can be stable on the model and oscillate on the bench.

Two engineering habits address this. The classical one is margins: keep 6 dB of gain margin and 45° of phase margin and hope the unmodelled effects are smaller. The modern one, introduced here, makes the hope quantitative. First, a single number, the **sensitivity peak Ms**, guarantees a gain *and* a phase margin simultaneously. Second, the model error is bounded by a frequency-dependent weight W, and the **small-gain test** ‖WT‖∞ < 1 then *guarantees* stability for every plant within that bound. Third, a theorem (the Bode sensitivity integral) tells you what no controller can do, so you stop trying.

Success criteria, written first: the H∞ norm by sweep and by Hamiltonian bisection agree to 10⁻⁶; every measured margin satisfies the bound its Ms implies; the weight covers the true actuator's relative error at every frequency; the designs that pass the small-gain test are stable on the true actuator, and for a design that fails, an explicit perturbation of size below 1 destabilises the loop; the Bode integral is reproduced to 10⁻⁴.

## 2. Understand and derive

### 2.1 The four closed-loop maps

Write the loop with controller C, plant G and loop transfer L = CG. Add a reference r, a disturbance d at the output and noise n on the measurement. Solving the loop equations y = G u + d, u = C(r − y − n):

y = T(r − n) + S d,  r − y = S(r − d) + T n,  u = KS (r − n − d),   (1)

with **S = 1/(1 + L)** (the sensitivity), **T = L/(1 + L)** (the complementary sensitivity) and KS = C/(1 + L) (control effort per unit reference); a fourth map GS = G/(1 + L) carries input disturbances to the output. Two facts organise everything:

- **S + T = 1** at every frequency. Disturbance rejection wants |S| small; noise rejection and robustness (2.4) want |T| small. Both cannot be small at the same frequency, so good loops make S small at low frequency and T small at high frequency, and cross over in between.
- For a type-1 loop (an integrator in L), S(0) = 0 and T(0) = 1: zero steady-state error to a step (`P09.08`).

As polynomials, with L = n_L/d_L and the characteristic polynomial χ = d_L + n_L: S = d_L/χ, T = n_L/χ, KS = n_C d_G/χ, GS = n_G d_C/χ. All four share the closed-loop poles (the roots of χ); the loop is internally stable when they are all in the left half plane (`is_closed_loop_stable`).

### 2.2 The H∞ norm, two ways

For a stable transfer function G, the **H∞ norm** is its peak gain, ‖G‖∞ = sup_ω |G(jω)| (for a matrix, the largest singular value). Ms = ‖S‖∞ and Mt = ‖T‖∞ are the two peaks of the loop.

**By a sweep.** Evaluate |G(jω)| on a logarithmic grid, take the largest, and refine with a golden-section search in log ω around it. This can miss a very narrow peak between grid points; for a lightly damped resonance the peak's width is about 2ζω_n.

**By the Hamiltonian.** For G = C(sI − A)⁻¹B with D = 0 and A Hurwitz, form

H(γ) = [[A, BBᵀ/γ²], [−CᵀC, −Aᵀ]]   (2)

(with a nonzero D the blocks change; the code uses the general form). Then **jω is an eigenvalue of H(γ) if and only if γ is a singular value of G(jω)** (Problem 2d). So ‖G‖∞ ≥ γ exactly when H(γ) has an imaginary-axis eigenvalue, and bisection on γ finds the norm to any tolerance without a frequency grid. This is the Boyd–Balakrishnan / Bruinsma–Steinbuch idea used in control toolboxes. For the lesson's chosen servo design both methods give Ms = 1.18243359 at 71.2 rad/s; on random systems they agree to 10⁻⁶ (the test), and for a resonance ω_n²/(s² + 2ζω_ns + ω_n²) both reproduce the closed form 1/(2ζ√(1 − ζ²)) to 10⁻⁷.

### 2.3 The disk margin from Ms

The distance from the Nyquist curve L(jω) to the critical point −1 is |1 + L(jω)| = 1/|S(jω)|, so

**min_ω |1 + L(jω)| = 1/Ms**:   (3)

the curve keeps out of a disk of radius 1/Ms centred at −1. That disk cuts the real axis between −1 − 1/Ms and −1 + 1/Ms. A gain change k moves real-axis crossings by the factor k, and cannot bring one onto −1 while k stays in

**(Ms/(Ms + 1), Ms/(Ms − 1)).**   (4)

At the gain crossover, |L| = 1 and L lies outside the disk; the unit circle meets the disk's edge where the angle from −1 is 2 arcsin(1/(2Ms)), so

**PM ≥ 2 arcsin(1/(2Ms)).**   (5)

These bounds hold *together*: a loop with Ms = 2 tolerates a gain between 2/3 and 2 (−3.5 dB to +6.0 dB) and 29.0° of phase, even if both are lost at once. That is why Ms between 1.2 and 2 is a common robustness specification. LQR's guarantee (`P09.16`) is the special case Ms ≤ 1 at the plant input: (½, ∞) and 60°.

### 2.4 Multiplicative uncertainty and the small-gain test

Write the real plant as **G = G₀(1 + W Δ)**, where W(s) is a known stable weight that bounds the relative error and Δ is any stable perturbation with ‖Δ‖∞ ≤ 1. The relative error of the true servo, Δ_true = G_true/G₀ − 1 = ω_a²/(s² + 2ζ_aω_as + ω_a²) e^{−0.005s} − 1, is near zero at low frequency and rises to a peak of 2.463 at 152.5 rad/s (the compliance resonance). The weight

W(s) = (s/40 + 0.05)/(s/160 + 1)

(5 % uncertainty at low frequency, 100 % near 40 rad/s, 400 % at high frequency) covers it: max_ω |Δ_true|/|W| = 0.905 ≤ 1.

With C stabilising G₀, the perturbed return difference factors as

1 + CG = 1 + L₀ + L₀WΔ = (1 + L₀)(1 + T W Δ).   (6)

The first factor has no right-half-plane zeros (the nominal loop is stable). The second cannot encircle the origin if |TWΔ| < 1 at every frequency, and |TWΔ| ≤ |WT| ‖Δ‖∞ < 1 when **‖WT‖∞ < 1**. So the small-gain test guarantees stability for *every* plant in the set: this is **robust stability**. The test is also necessary: if |W T(jω*)| ≥ 1 at some ω*, choose Δ(jω*) = −1/(WT(jω*)), whose size is ≤ 1, and realise it as a stable all-pass c(a − s)/(a + s) with the right phase (`destabilising_perturbation`). Then 1 + TWΔ = 0 at jω*: a closed-loop pole on the imaginary axis.

The test is honest about its question, "stable for every Δ in the set?", and that question is more demanding than "stable for the actual plant?". The gap between the two is the price of not knowing the phase of the model error.

### 2.5 The Bode sensitivity integral: what no controller can do

If L has relative degree at least two and the open loop has right-half-plane poles p_i, then for any stabilising controller

**∫₀^∞ ln |S(jω)| dω = π Σ Re(p_i).**   (7)

For a stable open loop the area of ln|S| is zero: every region where |S| < 1 (disturbances attenuated) is paid for by a region where |S| > 1 (amplified). Pushing S down in one band pushes it up elsewhere, like a waterbed. An unstable open loop pays more: πp extra. For the representative unstable airframe G = 20/((s − 2)(0.05s + 1)), divergence rate p = 2 s⁻¹, the integral is 2π = 6.2832 whatever stabilising controller is used. If the controller must keep |S| ≤ ε below ω₁ and the actuator cannot act above ω₂ (|S| ≈ 1 there), (7) forces

ln Ms ≥ (πp − ω₁ ln ε)/(ω₂ − ω₁),   (8)

a lower bound on the sensitivity peak that depends only on the physics (Problem 3).

### 2.6 Where the model stops being valid

- **Unstructured Δ.** The small-gain test allows any complex Δ; a real gain error is a very special Δ, and the test can be very conservative for it (Problem 5).
- **One weight, one actuator.** The weight was fitted to one unit. Production spread, temperature and wear move the compliance mode; the weight must cover all units.
- **Padé.** The delay is exact in frequency-domain evaluations and a Padé approximant only where a polynomial is needed (stability tests); the fourth-order approximant matches the exact relative error to 4 × 10⁻⁷ up to 300 rad/s and 8 × 10⁻⁴ at 1000 rad/s.
- **Linear.** Saturation and backlash are not uncertainty in this sense; they need `P09.19`.

## 3. Build the computational model

Open `starter/main.py`. The frequency response, the Padé approximant, the plant, true plant and weight models, the lead compensator and the design study are given.

### Step 1: the loop algebra

`closed_loop_tfs(nC, dC, nG, dG)` builds n_L = n_Cn_G, d_L = d_Cd_G, χ = d_L + n_L with `np.polymul` and `np.polyadd` and returns the five pairs. `is_closed_loop_stable` takes `np.roots` of χ. `tf_to_ss` writes the controllable canonical form: normalise the denominator, D = leading numerator coefficient, C from the remainder num − D·den.

### Step 2: the H∞ norm

```python
def hamiltonian_has_imaginary_eigs(A, B, C, D, gamma, tol=1e-8):
    R = gamma**2 * np.eye(D.shape[1]) - D.T @ D
    Rt = gamma**2 * np.eye(D.shape[0]) - D @ D.T
    F = A + B @ np.linalg.solve(R, D.T @ C)
    H = np.block([[F, B @ np.linalg.solve(R, B.T)],
                  [-gamma**2 * C.T @ np.linalg.solve(Rt, C), -F.T]])
    lam = np.linalg.eigvals(H)
    return np.any(np.abs(lam.real) < tol * max(1.0, np.abs(lam).max()))
```

`hinf_norm_bisection` refuses an unstable A, starts the lower bound at max(σ(D), σ(G(0))), doubles an upper bound until the test is false, and bisects to relative tolerance 10⁻⁸. `hinf_norm_sweep` evaluates σ_max(C(jωI − A)⁻¹B + D) on the grid, refines by 60 golden-section steps, and also checks ω = 0 and ω → ∞.

### Step 3: margins and uncertainty

`disk_margin_bounds(Ms)` is (4) and (5). `margins(w, L)` finds the first fall of |L| through 1 (phase margin, with the unwrapped phase) and every sign change of Im L with Re L < 0 (gain margins above 1 and below 1), interpolating linearly; it also counts a loop whose lowest-frequency sample is on the negative real axis, as for an unstable plant without an integrator. `robust_stability_index` returns ‖WT‖∞ by bisection. `destabilising_perturbation(w*, target)` returns c(a − s)/(a + s) with |c| = |target| and a chosen so the phase matches (the factor (a − jω)/(a + jω) sweeps the phase from 0 to −180°, and the sign of c supplies the other half-plane). `bode_sensitivity_integral` integrates ln|S| on a grid uniform in log ω over 12 decades below w_max = 10⁶ rad/s and adds the tail w_max ln|S(jw_max)|, exact when ln|S| decays as 1/ω².

### Verify

Run `python -m pytest curriculum/phase-09/lesson-18-robustness-uncertainty`. The tests:

1. **Realisation**: `tf_to_ss` against `scipy.signal.freqs` and the poles of the denominator; improper input refused.
2. **Loop identities**: S + T = 1, S = 1/(1 + L), KS = CS, GS = GS at 40 000 frequencies; S(0) = 0, T(0) = 1, KS(∞) = K_cp/z; stability verdicts for four loops.
3. **H∞, the independent check**: sweep against bisection on random systems of order 1, 2, 3 and 5 to 10⁻⁶; three resonances and a first-order lag against their closed forms; the Hamiltonian test true at 0.99‖G‖ and false at 1.01‖G‖; unstable A refused.
4. **Disk margin, the independent check**: GM = 3 for 2/(s(s + 1)(s + 2)); the formula at Ms = 2 and Ms → 1; on five loops (two nominal servos, two true servos with the exact delay, the unstable airframe), the measured PM, upper GM and lower GM satisfy the bounds from each loop's own Ms; the airframe's lower GM is exactly 1/6.
5. **Small gain**: the weight covers the true error; the Padé model equals the exact delay model to 10⁻⁶ below 30 rad/s; ‖WT‖∞ < 1 and 20 random all-pass and constant perturbations of size up to 0.99 are all stable; for ‖WT‖∞ > 1 the constructed Δ is all-pass of size below 1 and puts a closed-loop root within 10⁻⁶ω* of jω*; the design table orders ‖WT‖∞ and true stability with gain.
6. **Bode integral**: 0 to 10⁻³ for the servo, 2π to 10⁻⁴ for the airframe with two controllers.

### Compare

`scipy.signal` has no H∞ norm, no margin function and no sensitivity tools; the oracles here are `freqs` (for the realisation) and closed forms. Control toolboxes compute the H∞ norm with refinements of the same Hamiltonian idea (iterations that converge faster than bisection) and margins by the same crossing search. Our sweep and bisection agree to eight digits on Ms of the chosen design; the sweep is faster on this small problem and the bisection is grid-free, which is what you want when a lightly damped mode could hide between grid points.

Running `main()` prints the study:

```
servo G0 = 30/(s(0.04 s + 1)); unmodelled mode 150 rad/s, zeta 0.3, delay 5 ms: peak |Delta| = 2.463 at 152.5 rad/s
weight W = (s/40 + 0.05)/(s/160 + 1): max |Delta|/|W| = 0.905 (<= 1: covered)
chosen design Kc = 0.8: Ms by sweep 1.18243359 at 71.22 rad/s, by Hamiltonian bisection 1.18243359
   ||KS||inf = 3.200 V/rad (peak control effort per unit reference), ||T||inf = 1.0000
Kc = 2.5: ||WT|| = 1.594 at 93.2 rad/s; an all-pass Delta of size 0.627 puts a closed-loop pole at -0.0000+93.1671j
Bode integral: servo 1.80e-05 (theory 0); unstable airframe 6.2832 (theory pi p = 6.2832); its Ms = 1.995
```

with the design table of part 4, and writes `robustness.png` (|T| of five designs under 1/|W| with |Δ_true|; ln|S| of the servo and the airframe).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running.

1. **Design sweep.** For lead compensators (K_c, z, p) = (0.5, 20, 80), (0.8, 20, 80), (1.2, 20, 80), (1.6, 25, 100), (2.5, 25, 150) and (3.2, 25, 200), record the nominal crossover and PM, Ms, Mt, the disk bounds, ‖WT‖∞, and the true servo's stability, Ms, PM and GM (`design_table`). Predict: which designs pass the small-gain test, and whether any design that fails it is nevertheless stable on the true actuator.
2. **The conservatism.** For the (z, p) = (20, 80) family find the largest K_c that passes the small-gain test and the largest that is stable on the true servo (Problem 6). Predict the ratio.
3. **The waterbed.** For the unstable airframe with three stabilising leads, compute the Bode integral and Ms. Predict the integral before running.

Raw data format: CSV with header `Kc, z, p, wgc_rad_s, pm_nom_deg, Ms, Mt, gm_low_bound, gm_high_bound, pm_bound_deg, WT_inf, true_stable, Ms_true, pm_true_deg, gm_true`; `Kc_robust, Kc_true_limit`; `Kc, z, p, bode_integral, Ms`. Each header line names the plant, the true actuator and the weight.

What you should find:

| K_c, z, p | ω_gc (rad/s) | PM nom. (°) | Ms | ‖WT‖∞ | true stable | true Ms | true PM (°) | true GM | bounds from true Ms: GM ≥, PM ≥ |
|---|---|---|---|---|---|---|---|---|---|
| 0.5, 20, 80 | 15.86 | 84.8 | 1.121 | 0.414 | yes | 1.331 | 76.4 | 4.87 | 4.02, 44.1° |
| 0.8, 20, 80 | 25.95 | 78.3 | 1.182 | 0.658 | yes | 1.610 | 64.0 | 3.05 | 2.64, 36.2° |
| 1.2, 20, 80 | 38.34 | 70.0 | 1.257 | 0.974 | yes | 2.181 | 47.1 | 2.03 | 1.85, 26.5° |
| 1.6, 25, 100 | 43.94 | 66.3 | 1.263 | 1.108 | yes | 2.574 | 39.4 | 1.74 | 1.64, 22.4° |
| 2.5, 25, 150 | 68.26 | 65.5 | 1.272 | 1.594 | yes | 21.29 | 8.0 | 1.05 | 1.05, 2.7° |
| 3.2, 25, 200 | 87.89 | 66.3 | 1.263 | 1.887 | **no** | — | −58.5 | — | — |

The nominal columns are almost useless: every design has Ms ≈ 1.2–1.3 and a nominal phase margin above 65°, and the last one is unstable on the real actuator. The nominal model cannot see the compliance mode; the weight can. The small-gain test passes the first three designs and fails the last three; the true servo is stable for the first five. So the test was right for every design it passed, and cautious for two it failed: the K_c = 1.6 design is stable with 39° of true phase margin, and the K_c = 2.5 design is stable but with 8° and a true Ms of 21, which no engineer would fly. For that design the explicit perturbation has size 0.627 and puts a pole at 93.17j rad/s: a plant inside the uncertainty set that the controller cannot handle, even though the particular plant on the bench is (barely) stable. In every stable row the measured margins of the true loop satisfy the disk bounds from its own Ms, to the third digit in the 2.5 row (GM 1.05 against the bound 1.05). Problem 6 puts numbers on the conservatism: the robust gain limit is K_c = 1.234, the true stability limit 2.436, a factor 2. In the waterbed experiment the integral is 6.28319 for all three controllers (Ms 1.995, 1.947 and 2.018), against 2π = 6.28319.

## 5. Validate and improve

| quantity | theory (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| Ms, chosen design | ‖S‖∞ | sweep 1.18243359 | bisection 1.18243359 | < 10⁻⁸ | — |
| resonance peak, ζ = 0.25 | 1/(2ζ√(1 − ζ²)) = 2.06559 | 2.06559 | both methods | < 10⁻⁷ | tolerance |
| disk bound, K_c = 2.5 true loop | GM ≥ Ms/(Ms − 1) = 1.05 | measured 1.05 | — | third digit | bound nearly tight |
| weight coverage | ≤ 1 | 0.905 | — | — | weight chosen to fit |
| destabilising Δ, K_c = 2.5 | 1 + WTΔ = 0 at 93.2 rad/s | pole at 93.1671j | `np.roots` of χ | Re < 10⁻⁴ | — |
| Bode integral, servo | 0 | 1.8 × 10⁻⁵ | — | 1.8 × 10⁻⁵ | grid and tail |
| Bode integral, airframe | 2π = 6.28319 | 6.28319 | three controllers | < 4 × 10⁻⁶ | grid and tail |
| robust vs true gain limit | — | 1.234 vs 2.436 | — | factor 2.0 | first-order weight, complex Δ |

The tools agree with each other and with the theorems to the precision claimed. The validation question is whether the weight is right, and that cannot be answered from the model: it needs measured frequency responses of several actuators. What to improve: a weight shaped like the resonance (Extension 1) would let more gain pass the test; a structured analysis (Extension 3) would distinguish the real gain error from the dynamic error. What not to improve: the sweep resolution. Bisection already gives the norm to 10⁻⁸.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (3)–(6); `robustness.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) choosing a lead gain for the servo, with its nominal and true margins, ‖WT‖∞, and the hardware test that would confirm or refute the weight.

### Adversarial self-checks

- **What if the compliance mode drops to 120 rad/s?** The weight still covers it, barely (max |Δ|/|W| = 0.997, Problem 5). A further drop needs a new weight and a lower gain.
- **What if the gain is doubled?** From K_c = 0.8 to 1.6 the nominal Ms barely changes (1.18 to 1.26), but ‖WT‖∞ goes from 0.66 to 1.11 and the true phase margin from 64° to 39°: the danger is invisible on the nominal model.
- **Which assumption fails first?** The single weight for all actuators: production spread and temperature move the mode. Next, linearity: gear backlash creates a dead zone that the small-gain test does not describe.
- **What measurement could falsify the analysis?** A swept-sine test of the real actuator: if |G_meas/G₀ − 1| exceeds |W| anywhere, the guarantee is void. A closed-loop test at the chosen gain with a gain-injection point would measure the true margins directly.
- **Does the result survive a change of method?** Sweep and Hamiltonian bisection agree to 10⁻⁸; the true-plant stability verdict is the same with Padé orders 3 to 6 (the true limit is K_c = 2.43646 for each).
- **Which decision changes if the uncertainty is twice as large?** Doubling W doubles ‖WT‖∞: only the K_c = 0.5 design (0.83) would still pass. The design choice would move from 0.8 or 1.2 down to 0.5, halving the crossover.

## Check yourself

1. Define S and T and prove S + T = 1. Which one should be small at low frequency, and why?
2. What is the H∞ norm of a stable transfer function, and what is it for 3/(s + 2)?
3. A loop has Ms = 1.5. What gain range and phase margin are guaranteed?
4. State the small-gain test for multiplicative uncertainty and say what "robust stability" means.
5. Why is ‖WT‖∞ < 1 sufficient for robust stability, and why is it also necessary?
6. Give the Bode sensitivity integral for a stable open loop and for one with a right-half-plane pole p.
7. Why can a design with a nominal phase margin of 66° be unstable on the real actuator?
8. Why is the small-gain test conservative for a real gain error?

## References

- `@mit-ocw-16-323-principles-of-optimal-control` (supplementary-free): the robustness lectures (sensitivity functions, unstructured uncertainty, the small-gain theorem) in its lecture notes.
- `@scipy-userguide` (supplementary-free): `scipy.signal.freqs`, the oracle of the realisation test.
- The Hamiltonian bisection, the disk-margin geometry, the destabilising all-pass construction and the servo study are written here as original instruction.
