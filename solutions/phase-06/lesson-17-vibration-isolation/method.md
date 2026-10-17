# P06.17  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. k = 4 × 12 kN/m; `natural_frequency(k, m)`; `damping_coefficient`; `static_deflection(1, fn)`, `static_deflection(5, fn)`; check with √(g/δ)/(2π). 2. `transmissibility` at r = 1, 2, 500/fₙ; 1 − T. 3. `peak_transmissibility(0.1)` against 1/(2ζ).

## 2.
(a) Substitute exponentials, divide by k, take the modulus and the argument. For the force: the force through spring and damper to a fixed base is (k + icω)X with X = F₀/(k − mω² + icω). (b) With r² = 2 the numerator and denominator are both 1 + 8ζ²; for r² > 2, (1 − r²)² > 1, so the denominator exceeds the numerator. (c) d/du of (1 + 4ζ²u)/((1 − u)² + 4ζ²u) = 0 gives 2ζ²u² + u − 1 = 0; take the positive root; expand for small ζ. (d) Replace k + icω by k(1 + iη); for r ≫ 1, T_h ≈ √(1 + η²)/r² against T ≈ 2ζ/r.

## 3.
(a) `fn_for_attenuation(4 × 5.6, 0.15, 0.1)`; `stiffness_for_frequency`; `static_deflection(1, fn)`. (b) Solve T(r, 0.1) = 1.1 for r < 1 with `brentq`; fₙ ≥ f₁/r. Compare the two bounds. (c) T at f₁/fₙ and at 8f₁/fₙ with fₙ from (a). (d) Solve `transmissibility_hysteretic(f₄/fₙ, 0.2) = 0.15` for fₙ; compare. Think of devices that remove one frequency rather than all high frequencies.

## 4.
(a) `choose_isolator` with the finer grids and your sway limit; read off the four numbers; check which constraints are within a few per cent of their limits. (b) `synthesize_random(..., 32.5, seed=s)`, `simulate_sdof_base`, discard 0.5 s; RMS of the acceleration and g₀ × RMS of z; compare with `random_response`; normalise by u. (c) Compare grid steps and the binding constraints.

## 5.
(a) `random_response`, `peak_transmissibility`, `transmissibility(500/65, ζ)` for each ζ. (b) Build the breakpoint arrays; `psd_grms`; `random_response` at 65 Hz for each ζ. (c) Locate the energy of each spectrum relative to fₙ and √2 fₙ; reason with T at resonance and T ≈ 2ζ/r far above it.

## 6.
(a) Generate the data with the given seed; `curve_fit` each model with `sigma = 0.03 * T`, starting near 70 Hz; χ² = Σ((T − model)/σ)²; dof = 25 − 2. (b) `np.polyfit(log10 f, log10 T, 1)` on f > 300 Hz; predictions −1 and −2. (c) Compare χ²/dof; compute what the viscous model with ζ = η/2 predicts at 500 Hz against the hysteretic one.
