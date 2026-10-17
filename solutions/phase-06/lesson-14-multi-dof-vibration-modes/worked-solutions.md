# P06.14  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. All results are outputs of the model, not measurements.

## 1. A two-mass instrument stack (s = 15: m₁ = 275 kg)

**Assumptions.** Lumped masses, massless linear springs, fixed base, axial motion.

(a) M = diag(275, 300) kg; K = [[7 × 10⁸, −2 × 10⁸], [−2 × 10⁸, 2 × 10⁸]] N/m. The characteristic equation 82 500 ω⁴ − 2.65 × 10¹¹ ω² + 10¹⁷ = 0 has roots ω² = 4.3674 × 10⁵ and 2.7754 × 10⁶ s⁻²: **ω = 660.86 and 1665.95 rad/s, f = 105.180 and 265.144 Hz**.

(b) x₂/x₁ = (k₁ + k₂ − ω²m₁)/k₂: **+2.8995** (in phase) and **−0.31615** (out of phase). Mass-normalised: φ₁ = (0.018908, 0.054823), φ₂ = (0.057261, −0.018103). φ₁ᵀMφ₂ = **6.8 × 10⁻¹⁷** (round-off).

(c) Γ = ΦᵀM(1, 1) = **21.647 and 10.316**; effective masses **468.58 and 106.42 kg** (81.5 % and 18.5 %), summing to **575.00 kg** = m₁ + m₂.

(d) `modes` and `participation` return the same values.

**Check.** If k₁ → ∞ the platform is locked and the instrument alone has ω = √(k₂/m₂) = 816.5 rad/s, between the two frequencies, as the interlacing property of symmetric eigenproblems requires.

**Common mistakes.** Forgetting k₂ in K₁₁; normalising to unit length instead of unit modal mass.

## 2. The modal toolkit

(a) From Kφᵢ = ωᵢ²Mφᵢ, pre-multiplying by φⱼᵀ: φⱼᵀKφᵢ = ωᵢ²φⱼᵀMφᵢ. Swapping i and j: φᵢᵀKφⱼ = ωⱼ²φᵢᵀMφⱼ. K and M symmetric make the left sides equal and the M terms equal, so (ωᵢ² − ωⱼ²)φᵢᵀMφⱼ = 0: for distinct frequencies **φᵢᵀMφⱼ = 0**, and then φᵢᵀKφⱼ = 0. For a repeated frequency every combination of the two modes is a mode; one can always choose an M-orthogonal pair (the Cholesky route and `eigh` do so automatically).

(b) The modes form a basis, so r = Φc for some c; multiplying by ΦᵀM gives c = ΦᵀMr = Γ. Then **rᵀMr = cᵀΦᵀMΦc = ΓᵀΓ = Σ Γᵢ²**.

(c) For ÿ = A e^{iωt}, (10) gives η̂ᵢ = −ΓᵢA/(ωᵢ² − ω² + 2iζᵢωᵢω). The absolute acceleration is ẍ = A r − ω²Φη̂, and the force the interface supplies is the total inertia force F = rᵀMẍ = A rᵀMr − ω² rᵀMΦη̂ = A ΣΓᵢ² − ω² Σ Γᵢη̂ᵢ, using rᵀMΦ = Γᵀ and rᵀMr = ΣΓᵢ². Term by term:
F/A = Σᵢ Γᵢ² [1 + ω²/(ωᵢ² − ω² + 2iζᵢωᵢω)] = Σᵢ Γᵢ² (ωᵢ² + 2iζᵢωᵢω)/(ωᵢ² − ω² + 2iζᵢωᵢω), which with rᵢ = ω/ωᵢ is **Σ m_eff,i (1 + 2iζᵢrᵢ)/(1 − rᵢ² + 2iζᵢrᵢ)**.

(d) ΦᵀCΦ = aΦᵀMΦ + bΦᵀKΦ = aI + b diag(ωᵢ²); the i-th diagonal entry is 2ζᵢωᵢ, so **ζᵢ = a/(2ωᵢ) + bωᵢ/2**. ζ(ω) is a U-shaped curve, passing through the two fitted points, with a minimum between them at ω = √(a/b): modes between the fitted frequencies are slightly *less* damped than the targets, modes above (and below) are *more* damped, increasingly so.

**Common mistakes.** In (c), summing the modal accelerations without the rigid-body term A r; in (d), assuming the third mode gets the same ζ as the first two.

## 3. The stack after separation

(a) K_ff = [[4, −4, 0], [−4, 6.5, −2.5], [0, −2.5, 2.5]] × 10⁸ N/m. K_ff(1, 1, 1) = 0: the null vector is **(1, 1, 1)**, a **rigid-body translation** that stretches no spring.

(b) Frequencies **0, 160.05 and 295.74 Hz**; shapes (normalised to x₁ = 1): (1, 1, 1); (1, 0.6208, −1.4927); (1, −0.2949, 0.0769). The fixed-base stack had 73.54, 180.04 and 420.21 Hz. Removing k₁ changes the boundary condition, so every mode changes: the first fixed-base mode (the stack bouncing on its interface spring) becomes the rigid-body mode, and the free–free flexible modes are new shapes (masses moving against each other about the centre of mass), not the old ones with a mode deleted.

(c) Σmᵢφᵢ = **33.17** for the rigid mode (its mass-normalised Γ, √1100) and **7.9 × 10⁻¹⁵ and 2.6 × 10⁻¹⁵** for the flexible ones: zero momentum. Orthogonality with the rigid mode, φ_rigidᵀMφ = (1, 1, 1)Mφ = Σmᵢφᵢ = 0, *is* the statement that the flexible modes do not move the centre of mass.

(d) A thruster pulse excites the rigid-body mode (the stack accelerates: Δv) and every flexible mode whose shape has a non-zero entry at the thruster's location; the flexible modes then ring and decay with their damping, while the rigid-body mode keeps its new velocity (zero frequency: no restoring force).

**Common mistakes.** Using `np.linalg.cholesky` on the singular K; reading the free–free frequencies as the fixed-base ones shifted.

## 4. Thrust onset (s = 15: a₀ = 2.5 g = 24.52 m/s²)

(a) Each spring carries everything above it: **26 968, 23 291 and 8 581 N** (1100, 950 and 350 kg × a₀).

(b) Peak spring forces over 0.1 s (N) and dynamic amplification factors:

| model | k₁ (interface) | k₂ | k₃ |
|---|---|---|---|
| direct | 49 656 (1.841) | 45 602 (1.958) | 19 596 (2.284) |
| 1 mode | 47 921 (−3.5 %) | 45 363 (−0.5 %) | 20 602 (+5.1 %) |
| 2 modes | 48 968 (−1.4 %) | 46 053 (+1.0 %) | 19 505 (−0.5 %) |
| 3 modes | 49 656 (< 10⁻⁸) | 45 602 (< 10⁻⁸) | 19 596 (< 10⁻⁸) |

(c) Each mode's contribution peaks at up to twice its static share, and the static force in a spring is the sum of all modes' static shares. If two modes' contributions have opposite signs statically (one adds, one subtracts) but their oscillating parts add at some instant, the peak exceeds twice the total. Here the payload spring's static force is the small difference of large modal parts, and it reaches 2.28 times static.

(d) The payload's peak absolute acceleration is **5.73 g**, 2.29 times a₀: consistent with the k₃ spring force (the only force on the payload).

**Common mistakes.** Using only mode 1 because it carries 92 % of the mass: the interface load is then 3.5 % low and the payload spring 5 % high; for a step input, which excites every mode, keep them all.

## 5. A damper that is not proportional (s = 15: c_d = 4 × 10⁴ N·s/m)

(a) Diagonal modal ratios **(0.005882, 0.07577, 0.000776)**; exact from the state matrix **(0.005853, 0.07582, 0.000771)**, damped frequencies 73.573, 179.457, 420.187 Hz. The diagonal approximation is within 0.5 % in ζ for this light, well-separated case.

(b) Payload driving-point receptance peak near f₁: direct 6.0764 × 10⁻⁷ m/N, modal with diagonal ζ 6.0559 × 10⁻⁷ m/N, **ratio 0.9966**.

(c) 1 g at f₁: payload **105.5 g**, interface force **844 kN**: with only this damper, mode 1 has ζ = 0.0059 and the resonance is very sharp (a sweep from 60 to 90 Hz peaks at 105.8 g).

(d) A damper dissipates energy in proportion to the square of its stretch. In mode 1 the service and payload modules move nearly together (φ₃ − φ₂ = 0.0117 against 0.0655 in mode 2), so the damper hardly works in mode 1 and works hard in mode 2. Mode 1 stretches the adapter–service spring k₂ most (relative motion 0.0160 and the largest share of the mode's strain energy): a damper across k₂ would give mode 1 about (0.0160/0.0117)² ≈ 1.9 times more damping for the same c_d. Put dampers where the mode you want to damp stores its strain energy.

**Common mistakes.** Reading the coupling measure as a percentage error; using undamped modes for the state-space eigenvalues.

## 6. Modal models in a test campaign

(a) For the notch at f₁, either model gives the frequency; what sets the notch level is the response at resonance, which both models compute from an *assumed* ζ. Use the more detailed model for mode shapes and effective masses (it has the local modes and the lateral coupling), but set the notch from a measured low-level sweep scaled to full level, not from either model.

(b) Testable: (1) the damping is higher than assumed (rerun with ζ from the low-level sweep: the force scales as 1/ζ); (2) the interface is not rigid (add the shaker or fixture compliance as a fourth spring and see whether the predicted force falls). Not testable by the model alone: the test article differs from the model (a missing mass, a loose joint): only measured mass properties and a modal survey settle it.

(c) The launch vehicle's coupled-loads analysis needs a reduced model of the spacecraft at its interface, and the fundamental frequencies and effective masses are what decide whether the two structures couple and how large the interface loads are; they also hide proprietary design detail. (In practice a reduced "Craig–Bampton" model is exchanged.)

(d) The modal model is a prediction; its frequencies may be a few per cent off and its damping is a guess. A low-level sweep (a fraction of a g) measures the actual resonances and Q's safely, lets the notch be placed where the peak really is, and checks the model before the full-level run can overstress the hardware. It costs minutes; skipping it risks the hardware.

## Check-yourself answers

1. Because M multiplies ω² on the right-hand side; with M = LLᵀ and φ = L⁻ᵀu it becomes L⁻¹KL⁻ᵀu = ω²u, with a symmetric matrix.
2. φᵢᵀMφⱼ = 0 and φᵢᵀKφⱼ = 0 for ωᵢ ≠ ωⱼ; they follow from the symmetry of M and K.
3. It scales each mode so φᵢᵀMφᵢ = 1; then ΦᵀMΦ = I and ΦᵀKΦ = diag(ωᵢ²).
4. η̈ᵢ + 2ζᵢωᵢη̇ᵢ + ωᵢ²ηᵢ = φᵢᵀf(t), with ηᵢ(0) = φᵢᵀMx(0) and η̇ᵢ(0) = φᵢᵀMẋ(0).
5. ΦᵀCΦ must be diagonal; C = aM + bK guarantees it. A mode far above the fitted ones gets ζ ≈ bω/2, increasingly large (stiffness-proportional damping dominates).
6. Γᵢ = φᵢᵀMr measures how strongly a rigid base motion r drives mode i; the effective mass is Γᵢ² (mass-normalised modes), and the effective masses sum to rᵀMr, the total mass for translation.
7. Well below their own frequencies, the omitted modes respond quasi-statically: their masses move with the base, so they contribute their effective mass times the base acceleration, which is exactly a rigid residual mass.
8. About m_eff,i/(2ζᵢ) times the base acceleration, i.e. (m_eff,i/m_total)/(2ζᵢ) times the rigid-body force: 23 times for the lesson's stack at ζ = 0.02.
