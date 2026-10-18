# P10.10  Prandtl Lifting-Line Theory

> Each section of a wing behaves like an airfoil, but at the angle the wake lets it see. Write the circulation as a Fourier series in the span angle and Prandtl's integral equation becomes a small linear system.

| | |
|---|---|
| **Learning outcomes** | 1. Derive Prandtl's lifting-line equation from the section lift law and the downwash of the trailing vortex sheet, and Glauert's monoplane equation from the substitution y = −(b/2) cos θ and Γ = 2bV∞ Σ Aₙ sin nθ. 2. Derive CL = πAR A1, CDi = πAR Σ nAₙ², e = 1/(1 + δ), and Glauert's lift-slope factor τ. 3. Implement and test a collocation solver for any chord, twist, section lift slope and zero-lift angle, symmetric or not, and verify it against the elliptic closed form CL = a0α/(1 + a0/(πAR)), e = 1. 4. Measure the observed order of the collocation convergence and explain why it is 4 for the rectangular wing and 2 for a tapered one. 5. Report the trainer's lift-curve slope, span efficiency and induced drag at cruise, and how δ and τ depend on taper and aspect ratio. |
| **Prerequisites** | `P10.09` (horseshoes, downwash, induced drag, the elliptic load), `P10.04` (section lift slope a0 and α_L0), `P02.09` (Fourier series), `P01.13` (linear systems), `P03.11` (observed order). Skill check: you can show ∫₀^π sin nθ sin mθ dθ = (π/2)δ_nm. |
| **Engineering question** | What lift-curve slope, span efficiency and induced drag does the reference trainer's rectangular wing (b = 10 m, c = 1.2 m, NACA 2412) have at its 55 m/s cruise? Answer, by lifting line with the thin-airfoil section (a0 = 2π, α_L0 = −2.077°): CL_α = 4.878 per radian (0.0851 per degree), 3.7 % below an elliptic wing's 5.067 (τ = 0.200); cruise CL = 0.331 at α = 1.81°; e = 0.934 (δ = 0.0707); CDi = 0.00448, an induced drag of 99.6 N, 7 % more than an elliptic load would give. |
| **Assumptions and limits** | Lifting-line theory: a large aspect ratio, an unswept straight lifting line, a flat wake, small induced angles, each section acting as a 2-D airfoil at its effective angle with a linear lift curve (no stall), steady incompressible inviscid flow. It fails for low aspect ratio (below about 4), sweep, and near the stall; it ignores the chordwise extent of the tips (`P10.11` measures that difference). a0 = 2π (thin airfoil) unless stated; real sections have a slightly lower slope. |
| **Analytical work** | The lifting-line integro-differential equation; Glauert's substitution and the evaluation of the downwash integral term by term; the monoplane equation; CL, CDi, δ and τ; the elliptic case solved exactly; the rolling-moment coefficient Cl = (πAR/4)A2; the regularity argument for the convergence order. |
| **Computational deliverable** | `starter/main.py`: `collocation_points`, `solve_monoplane`, `wing_coefficients`, `span_loading`, `rolling_moment_coefficient`, `lift_slope_and_tau`, `elliptic_closed_form`, `alpha_for_cl`, `convergence_study` (`Wing`, `rectangular_wing`, `tapered_wing`, `elliptic_wing`, `linear_twist` and `main` are given). Tests prove: the elliptic wing is exact (A_{n>1} < 10⁻¹², e = 1, CL and CDi to 10⁻⁶) for three AR, a0, α; linearity in α and zero lift at α_L0; a symmetric wing has no even coefficients and a symmetric load that vanishes at the tips; CL and CDi equal the span integrals of Γ and Γα_i; convergence at order 4 (rectangular) and 2 (tapered); δ has its minimum at taper 0.3–0.45, and the trainer's rectangular wing has δ = 0.0707 from an independent solve written in the test; CL_α → a0 as AR → ∞, τ = 0 for the elliptic wing; twist superposes and washout unloads the tips; anti-symmetric twist rolls the wing without changing CL, and Cl = (πAR/4)A2 equals the span integral of the load's moment. |
| **Practical activity** | **computational experiment**: a collocation convergence study with observed orders for three planforms, the elliptic closed-form check, a table of δ and τ against taper and aspect ratio, and the trainer's cruise analysis with a sensitivity to the section lift slope. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `lifting_line.png`, the derivation of the monoplane equation and of CL, CDi and δ, the convergence table, the δ–τ table, and a half-page note on the trainer's wing. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the monoplane equation is derived; all public tests pass; the observed orders are reported and explained; the note states CL_α, e and CDi with the assumptions behind them. |
| **References** | `@mit-ocw-16-01-unified-engineering-fluids` Fluids (spring) Lecture 6 "3D Vortex Filaments, Lifting Line Theory", Lecture 7 "Elliptical Lift Distribution" (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f07_sp/) and Lecture 9 "General Wings" (notes: General circulation distribution and downwash; Load distributions on typical planforms; Lift slope reduction) (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f09_sp/) (required-free). `@mit-ocw-16-100-aerodynamics` Lecture 17 "Prandtl's Lifting Line Introduction" and Lecture 18 "Force Calculations for Lifting Line" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre17_cg/) (required-free). `@leishman-intro-aerospace-flight-vehicles` "Lifting Line & Finite Wing Theory" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/lifting-line-theory/): Helmholtz's Vortex Theorems; Lifting Line Theory; Numerical Solution to the Monoplane Wing Equation; Wing Shape & Loading Distributions (required-free). `@ntrs-naca-tr-116-prandtl-modern-hydrodynamics` Part II: D. "Theory of the Monoplane" (https://ntrs.nasa.gov/citations/19930091180) (supplementary-free). `@ntrs-naca-tn-1422-sivells-three-wings` (Sivells, https://ntrs.nasa.gov/citations/19930082118) (required-free): measured against calculated characteristics of three aspect-ratio-9 wings, the comparison of extension 2. Optional textbooks for deeper study (optional-paid, not required; library or used copy): Anderson, *Fundamentals of Aerodynamics*, for the monoplane equation worked by hand; Katz & Plotkin, *Low-Speed Aerodynamics*, for lifting-line theory beside lifting-surface methods. |
| **Extensions** | ★ Replace the linear section lift law by a tabulated cl(α) with a stall and solve the nonlinear lifting line by fixed-point iteration (`P10.15` builds on this). ★ Compute the 3-D lift slope of the TN-1422 wing (AR 9, taper 0.4, 2° washout) and compare with the measured value in `curriculum/phase-10/data/finite_wing_cla.csv`, stating the a0 you used. ★ Write a discrete version from `P10.09`'s horseshoes (a "Weissinger-free" lifting line: horseshoes on the quarter-chord line, control points on it too) and compare with the Fourier solution. |
| **Time** | unmeasured |

## 1. The engineering problem

`P10.09` showed what a span load costs in induced drag once the load is known. A designer knows the wing (its chord, twist and sections) and wants the load. The two are coupled: the load sheds vorticity, the vorticity makes downwash, the downwash lowers each section's angle, and the angle sets the load. Prandtl's lifting-line theory (1918) closes that loop with one equation per spanwise station, and Glauert's Fourier formulation turns it into a small, well-conditioned linear system. It is still the quickest reliable tool for the span load of a straight wing, the source of the e and τ factors in every aircraft-design book, and the base of the nonlinear stall models of `P10.15` and the aeroelastic loads of later phases.

For the reference trainer this lesson answers concrete questions: how steep is the wing's lift curve (it sets the gust response and the static margin of later phases); what is its span efficiency (the induced drag at cruise and in the climb); where along the span is the local cl highest (where it will stall first; `P10.12`)? Success criteria, written first: the solver reproduces the elliptic closed form to 10⁻⁶ for any aspect ratio, section slope and angle; CL and CDi computed from the Fourier coefficients equal the span integrals of the load; the collocation error falls at a measured and explained order; the trainer's CL_α, e and CDi are reported with their assumptions and a sensitivity to the section slope.

## 2. Understand and derive

### 2.1 The lifting-line equation

Model the wing as a bound vortex Γ(y) along the quarter-chord line plus the trailing sheet of strength −dΓ/dy it sheds (`P10.09`). The downwash at the line is

w(y) = (1/4π) ∫_{−b/2}^{b/2} (dΓ/dη) dη/(y − η).   (1)

Each section is a 2-D airfoil (`P10.04`) at its **effective angle** α_eff = α + ε(y) − α_L0(y) − w/V∞, where ε is the geometric twist. Its lift per unit span is ρV∞Γ = ½ρV∞² c a0 α_eff, so

**Γ(y) = ½ V∞ c(y) a0(y) [α + ε(y) − α_L0(y) − w(y)/V∞].**   (2)

(1) and (2) together are Prandtl's integro-differential equation for Γ.

### 2.2 Glauert's Fourier solution

Substitute y = −(b/2) cos θ (θ = 0 at the left tip, π at the right) and

Γ(θ) = 2bV∞ Σ Aₙ sin nθ.   (3)

The sine series vanishes at both tips, as Γ must. The downwash integral (1) is evaluated term by term with Glauert's integral ∫₀^π cos nφ dφ/(cos φ − cos θ) = π sin nθ/sin θ:

w(θ)/V∞ = Σ nAₙ sin nθ / sin θ.   (4)

Putting (3) and (4) in (2), with μ(θ) = a0 c/(4b), gives the **monoplane equation**

**Σ Aₙ sin nθ (μn + sin θ) = μ (α + ε − α_L0) sin θ.**   (5)

Keep N terms and enforce (5) at N collocation stations θ_i: N linear equations for N unknowns, with no integral left to evaluate. For a symmetric wing the even terms vanish (they are anti-symmetric in y), so we keep n = 1, 3, 5, … and collocate on one half, θ_i = iπ/(2N), i = 1…N, the last at the root.

### 2.3 Forces from the coefficients

Lift: L = ρV∞∫Γ dy = ρV∞(b/2)∫₀^π Γ sin θ dθ; only n = 1 survives the orthogonality integral, so

**CL = πAR A1.**   (6)

Induced drag D_i = ρ∫Γw dy (`P10.09`) with (3) and (4): orthogonality again leaves only the squares,

**CDi = πAR Σ nAₙ² = CL²(1 + δ)/(πAR),  δ = Σ_{n≥3} n(Aₙ/A1)²,  e = 1/(1 + δ).**   (7)

δ ≥ 0: any departure from the elliptic load (A1 alone) costs induced drag. The rolling moment of an unsymmetric load comes from n = 2 alone: **Cl = (πAR/4)A2** (positive right wing down).

### 2.4 The elliptic wing (worked example)

If c(y) = c0√(1 − (2y/b)²) = c0 sin θ and α, α_L0 are constant, then μ = μ0 sin θ and (5) is satisfied by A1 alone: A1(μ0 + 1) sin θ = μ0 α_a sin θ with α_a = α − α_L0. With μ0 = a0c0/(4b) and AR = 4b/(πc0), μ0 = a0/(πAR), so A1 = α_a a0/(πAR + a0) and

**CL = a0 α_a /(1 + a0/(πAR)),  e = 1.**   (8)

For a0 = 2π: CL = 2πα_a AR/(AR + 2). For AR = 8.333 and α_a = 5°, CL = 0.44219. The solver must reproduce this with all higher coefficients at round-off.

### 2.5 The lift slope and τ

For a general flat untwisted wing, CL is still linear in α; Glauert wrote its slope in the elliptic form with a correction τ:

**CL_α = a0 / (1 + (a0/(πAR))(1 + τ)).**   (9)

τ (like δ) depends on planform and aspect ratio; for the elliptic wing τ = δ = 0. The code computes CL_α by solving at two angles and then reads τ from (9).

### 2.6 Why the convergence order depends on the planform

The collocation solution is a truncated Fourier series; its error is set by how fast the true coefficients Aₙ decay, which is set by how smooth Γ(θ) is as a *periodic* function of θ. A rectangular wing's chord, written in θ, is smooth (constant) and Γ(θ) is smooth across the root; its coefficients decay algebraically fast and the observed order is 4. A straight-tapered wing has c ∝ 1 − (1 − λ)|cos θ|, whose derivative jumps at the root (θ = π/2); one jump in the first derivative gives Aₙ ~ 1/n² behaviour and order 2. Any twist that is linear in |y| does the same. A convergence study that used only one planform would misreport the method.

### 2.7 Where the model stops

Lifting-line theory places all the lift on one line. That is exact in the limit AR → ∞ and acceptable for AR above about 6 on unswept wings, but it ignores the fact that the chord has length: near the tips the load falls over a distance comparable with the chord, which a lifting surface sees and a lifting line does not. `P10.11` computes the trainer's wing with a vortex lattice and finds a lift slope about 5 % lower; the difference falls like 1/AR. Sweep (the bound vortex is no longer perpendicular to the flow) and stall (the section law is no longer linear) are outside this lesson.

## 3. Build the computational model

Open `starter/main.py`. Given: the `Wing` dataclass (b, chord(y), twist(y), α_L0(y), a0; `area()`, `aspect_ratio()`), `rectangular_wing`, `tapered_wing`, `elliptic_wing`, `linear_twist`, `main`.

### Step 1: collocation and solve

`collocation_points(N, symmetric)` and `solve_monoplane(wing, α, N, symmetric)` assemble (5):

```python
n = 2 * np.arange(N) + 1 if symmetric else np.arange(1, N + 1)
th = collocation_points(N, symmetric)
y = -0.5 * wing.b * np.cos(th)
mu = wing.a0_at(y) * wing.chord(y) / (4 * wing.b)
M = np.sin(np.outer(th, n)) * (mu[:, None] * n[None, :] + np.sin(th)[:, None])
A = np.linalg.solve(M, mu * (alpha + wing.twist(y) - wing.alpha_L0(y)) * np.sin(th))
```

### Step 2: results

`wing_coefficients` is (6) and (7); `span_loading` returns Γ, the local cl = 2Γ/(V∞c) and α_i from (4); `rolling_moment_coefficient` is section 2.3; `lift_slope_and_tau` is (9); `elliptic_closed_form` is (8); `alpha_for_cl` uses linearity; `convergence_study` compares with a 1280-term reference and returns observed orders log(e_k/e_{k+1})/log 2.

Running `main()` prints (abridged):

```
elliptic wing AR 8.3333, alpha 5 deg: CL 0.44218658 (closed form 0.44218658), e 1.0000000000, max |A_n>1| = 9.2e-19
reference trainer (rectangular, b 10.0 m, c 1.2 m, AR 8.333, a0 = 2 pi, alpha_L0 = -2.077 deg):
  CL_alpha = 4.8779 /rad = 0.08513 /deg; elliptic would give 5.0671; tau = 0.2004
  cruise CL 0.3309 at alpha = 1.810 deg; e = 0.9339 (delta 0.0707); CDi = 0.004479 (elliptic 0.004183); induced drag 99.6 N
  y (m)    Gamma (m^2/s)   cl/CL    alpha_i (deg)
   0.0       12.262       1.1229    0.499
   2.5       11.735       1.0746    0.644
   4.0       10.044       0.9198    1.112
   4.5        8.344       0.7641    1.581
   4.9        4.555       0.4171    2.628
collocation convergence (alpha = 5 deg, flat wings, reference N = 1280):
  rectangular AR 8.33: |CL error| [1.59e-04 9.19e-06 5.56e-07 3.44e-08 2.14e-09 1.34e-10]
      observed orders CL [4.11 4.05 4.02 4.01 4.  ], CDi [4.14 4.07 4.03 4.01 4.  ]
  taper 0.4, AR 8.33: |CL error| [1.99e-03 5.26e-04 1.33e-04 3.32e-05 8.28e-06 2.05e-06]
      observed orders CL [1.92 1.99 2.   2.   2.02], CDi [1.77 1.96 2.   2.   2.02]
```

and writes `lifting_line.png`: cl/CL along the span for three planforms at the cruise CL, the convergence plot, and δ against taper.

### Verify

`python -m pytest curriculum/phase-10/lesson-10-prandtl-lifting-line` runs 13 tests (one parametrised over three cases), listed in the metadata. Two are independent of the coefficient formulas: CL and CDi recomputed by quadrature of Γ and of Γα_i along the span (the Kutta–Joukowski and near-field drag of `P10.09`) agree with πAR A1 and πAR Σ nAₙ² to 10⁻¹⁰ and 10⁻⁸.

### Compare

The elliptic closed form (8) is reproduced to eight digits, with higher coefficients at 10⁻¹⁸. Against the horseshoe stack of `P10.09`, which shares no code with this solver, the two-term load sin θ + a₃ sin 3θ gave e = 1/(1 + 3a₃²) to 2 × 10⁻⁴, the δ of (7). The established tool for this problem is a vortex-lattice code such as AVL; `P10.11` builds one and compares with this lesson's results.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Convergence.** CL error against a 1280-term reference for N = 5–160, three planforms. Predict each order from section 2.6. *Finding:* rectangular 4.11 → 4.00; taper 0.4: 1.92 → 2.02; elliptic chord with linear washout −4°: 2.00 (the twist's kink at the root sets it). The rectangular wing reaches 10⁻⁸ in CL with N = 40.
2. **δ and τ against taper and AR** (flat, untwisted, a0 = 2π):

| AR | | λ = 0 | 0.2 | 0.3 | 0.35 | 0.4 | 0.5 | 0.7 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|
| 4 | δ | 0.0949 | 0.0114 | 0.0057 | 0.0048 | 0.0048 | 0.0064 | 0.0135 | 0.0285 |
| | τ | 0.1777 | 0.0361 | 0.0256 | 0.0255 | 0.0277 | 0.0368 | 0.0658 | 0.1195 |
| 8.33 | δ | 0.1595 | 0.0261 | 0.0148 | 0.0133 | 0.0137 | 0.0181 | 0.0355 | 0.0707 |
| | τ | 0.2377 | 0.0565 | 0.0433 | 0.0444 | 0.0488 | 0.0650 | 0.1134 | 0.2004 |
| 12 | δ | 0.1938 | 0.0361 | 0.0220 | 0.0204 | 0.0213 | 0.0281 | 0.0533 | 0.1031 |
| | τ | 0.2651 | 0.0679 | 0.0545 | 0.0567 | 0.0628 | 0.0836 | 0.1442 | 0.2511 |

(AR 6: δ from 0.1293 to 0.0087 at λ = 0.35 to 0.0483 rectangular.) Predict where δ is smallest. *Finding:* at λ ≈ 0.35 for every AR (0.36 at AR 8.33, δ = 0.0133, e = 0.987); both factors grow with AR at fixed planform, because the higher the aspect ratio, the more the 2-D section slope makes the load follow the chord rather than the elliptic shape.
3. **The trainer at cruise.** *Finding:* CL_α = 4.878 per radian, α = 1.81° for CL = 0.331, e = 0.934, CDi = 0.00448 (99.6 N against 93 N elliptic). The local cl is 1.123 CL at the root and falls to 0.417 CL at y = 4.9 m: the rectangular wing is root-loaded in cl, so it stalls at the root first, a benign stall.
4. **Sensitivity to the section slope.** a0 = 0.9 × 2π (a viscous section is typically a few per cent below 2π; `P10.06`): CL_α = 4.480 per radian, e = 0.927. A 10 % lower a0 lowers the wing slope by 8 % and e by 0.7 points.
5. **An unsymmetric case.** ±2° anti-symmetric twist outboard of 3 m (an aileron-like input) on the flat wing at α = 5°: A2 = −0.00311, Cl = −0.0204 (right wing up), CL unchanged at 0.42567.

Raw data format: CSV tables `planform, N, err_CL, err_CDi` (1); `AR, taper, delta, tau` (2); `y, Gamma, cl_over_CL, alpha_i` (3); `a0_factor, CL_alpha, e` (4); plus `lifting_line.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| elliptic CL, AR 8.33, 5° | (8): 0.44218658 | 0.44218658 | < 10⁻⁸ | — |
| elliptic e | 1 | 1.0000000000 | < 10⁻¹⁰ | — |
| rectangular CL, N = 40 | reference N = 1280 | error 3.4 × 10⁻⁸ | — | truncation, order 4 |
| CL from ∫Γ dy | πAR A1 | equal | < 10⁻¹⁰ | — |
| CDi from ∫Γα_i dy | πAR ΣnAₙ² | equal | < 10⁻⁸ | — |

The solver is verified; its numerical error is negligible at N = 40. The model error is not: section 2.7's chordwise-tip effect (about 5 % on the trainer's lift slope, `P10.11`), the real section's lift slope (part 4, ±8 %), and the fuselage, which carries lift across the root and is absent here. The engineering answer: **the trainer's wing has CL_α = 4.88 per radian (0.085 per degree), e = 0.934 and, at the 55 m/s cruise, CDi = 0.0045 (100 N) at α = 1.8°; the induced drag is 7 % above an elliptic load's, and the local cl peaks at the root.** What to improve: use a measured or `P10.06` section slope; compare with the lattice of `P10.11`; for the TN-1422 wing, compare with measurement (extension 2).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (4)–(9); `lifting_line.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the trainer's CL_α, e and CDi with their assumptions.

### Adversarial self-checks

- **What if the aspect ratio doubles?** At fixed planform shape, δ and τ grow (rectangular: δ 0.0707 at AR 8.33, 0.1031 at 12), but CDi = CL²(1 + δ)/(πAR) still falls almost in proportion to 1/AR.
- **What if the angle goes to the stall?** The linear section law in (2) is then wrong; the theory keeps adding lift. Stall needs the nonlinear section curve (extension 1, `P10.15`).
- **Which assumption fails first?** For the trainer, the lifting line itself (a 5 % lift-slope difference from the lattice) and the section slope a0; for a swept or low-AR wing, everything.
- **What measurement could falsify the model?** A measured span load (pressure taps at several stations, as in TN-1422) that differs in shape, or a measured lift slope outside the a0 uncertainty.
- **Does the result survive a change of resolution?** Yes: N = 40 is converged to 3 × 10⁻⁸ in CL for the rectangular wing and to 3 × 10⁻⁵ for the tapered one.
- **Which decision changes if the uncertainty doubles?** If the lift slope were uncertain by ±16 % rather than ±8 %, the trim angle at cruise (1.8°) would be known only to about ±0.6° (against ±0.3°), enough to matter for a wing-incidence decision but not for performance.

## Check yourself

1. What are the two physical ingredients of Prandtl's lifting-line equation?
2. Why does the Fourier series Σ Aₙ sin nθ automatically satisfy the tip conditions?
3. Which coefficient sets the lift, and which the rolling moment?
4. Define δ and e. Why is δ never negative?
5. For a flat elliptic wing with a0 = 2π and AR = 6, what is CL_α?
6. Why does the rectangular wing's collocation solution converge at order 4 but a tapered wing's at order 2?
7. Where does a rectangular untwisted wing stall first, and why?
8. Name two situations where lifting-line theory should not be used.
