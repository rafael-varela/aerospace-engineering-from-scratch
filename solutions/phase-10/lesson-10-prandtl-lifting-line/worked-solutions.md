# P10.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown.

## 1. An elliptic wing by hand (s = 15: AR = 9)

**Assumptions.** Flat elliptic planform, a0 = 2π, linear section lift.

(a) CL_α = 2π/(1 + 2/9) = **5.141 /rad**; at 4° (0.06981 rad) **CL = 0.3589**; α_i = CL/(πAR) = 0.01269 rad = **0.727°**; CDi = CL α_i = **0.004556**.

(b) The solver gives CL = 0.35889, CDi = 0.0045556, A3 at round-off and **δ = 1.7 × 10⁻³¹**.

(c) Doubling the span at the same area quadruples AR, so CDi = CL²/(πAR) falls to **one quarter**.

**Common mistakes.** Using a0 in degrees; computing CDi with α instead of α_i.

## 2. Derivations

(a) With y = −(b/2) cos θ, dy = (b/2) sin θ dθ and dΓ/dφ = 2bV∞ Σ nAₙ cos nφ; (1) becomes w(θ) = (1/4π)∫₀^π 2bV∞ Σ nAₙ cos nφ dφ/((b/2)(cos φ − cos θ)) = (V∞/π) Σ nAₙ π sin nθ/sin θ = **V∞ Σ nAₙ sin nθ/sin θ**. In (2): 2bV∞ Σ Aₙ sin nθ = ½V∞ c a0 [α + ε − α_L0 − Σ nAₙ sin nθ/sin θ]. Divide by 2bV∞, call μ = a0c/(4b), multiply by sin θ: **Σ Aₙ sin nθ (μn + sin θ) = μ(α + ε − α_L0) sin θ**.

(b) L = ρV∞∫Γ dy = ρV∞²b² Σ Aₙ ∫₀^π sin nθ sin θ dθ = ρV∞²b²(π/2)A1, and CL = L/(½ρV∞²S) = πb²A1/S = **πAR A1**. D_i = ρ∫Γw dy = ρV∞²b² Σ_n Σ_m Aₙ mA_m ∫ sin nθ sin mθ dθ = ρV∞²b²(π/2) Σ nAₙ², so **CDi = πAR Σ nAₙ²**.

(c) c = c0 sin θ ⇒ μ = μ0 sin θ, and with A1 alone the left side is A1 sin θ (μ0 sin θ + sin θ) = A1 sin²θ (μ0 + 1), the right side μ0 α_a sin²θ: equal for every θ with **A1 = μ0 α_a/(1 + μ0)**. Since S = πbc0/4, AR = 4b/(πc0) and μ0 = a0c0/(4b) = a0/(πAR). CL = πAR A1 = πAR μ0α_a/(1 + μ0) = **a0α_a/(1 + a0/(πAR))**; δ = 0, **e = 1**.

(d) Rolling moment (positive right wing down) = −∫ y L′ dy = −ρV∞∫Γy dy. With y = −(b/2) cos θ: −ρV∞ × 2bV∞ × (−(b/2))(b/2) Σ Aₙ ∫ sin nθ cos θ sin θ dθ = (ρV∞²b³/2) Σ Aₙ ∫ sin nθ (½ sin 2θ) dθ = (ρV∞²b³/2)(π/4)A2. Divided by ½ρV∞²Sb: **Cl = (πAR/4)A2**.

**Common mistakes.** Losing the principal value in (1); forgetting that y increases with θ (left tip at θ = 0).

## 3. A partial-span flap

**Assumptions.** Flat rectangular wing, flapped sections modelled by α_L0 = −0.609 × 10° = −6.09° (thin-airfoil τ of a 25 % flap), α = 0.

(a) CL = ΔCL = **0.3442, 0.3400, 0.3421** for N = 40, 80, 160, and e (of this flap-only load) 0.666, 0.657, 0.662. The convergence is **oscillatory and slow**: the right-hand side of (5) jumps at the flap edges, and a Fourier series converges slowly and non-monotonically at a jump (the Gibbs phenomenon). Collocation points falling just inside or outside the edge change the answer at the 1 % level. ΔCL ≈ **0.34 ± 0.002**.

(b) Section increment 2πτδ = 0.668. Flapped-span rule: 0.6 × 0.668 = **0.401** (18 % high). With the 3-D slope ratio 4.878/2π: **0.311** (9 % low). The lifting line sits between the two quick rules: the flap's lift spreads beyond its edges (the induced flow of its edge vortices), which neither rule captures.

(c) cl = **0.357 at y = 2.9 m** (inside) and **0.234 at 3.1 m** (outside), with 0.547 at the root and 0.056 at 4.5 m. The model keeps Γ continuous but steep at the edge; a real wing sheds a concentrated **flap-edge vortex** there, with a local suction peak and noise, and the flap's side edge has 3-D separated flow. Lifting-line theory has no chordwise or side-edge detail.

**Common mistakes.** Taking the N = 40 value as converged; using the full span for the flapped-area rule.

## 4. Your own convergence study (s = 15: λ = 0.375)

(a) Errors 5.49 × 10⁻⁴, 1.38 × 10⁻⁴, 3.46 × 10⁻⁵, 8.63 × 10⁻⁶ for N = 10, 20, 40, 80; observed orders **1.99, 2.00, 2.00**: order 2, as the kink of c(y) at the root predicts.

(b) δ = **0.01332**, e = **0.9869**: λ = 0.375 is right at the minimum of δ (0.0133 at λ ≈ 0.36).

(c) **Order 2**: a linear twist ε ∝ |y| has a slope discontinuity at the root, which makes the right-hand side of (5) a kinked function of θ, exactly like a tapered chord (the lesson's elliptic-plus-washout case shows 2.00).

**Common mistakes.** Measuring the order against a reference that is not much finer than the finest N.

## 5. The lift slope of a measured wing

**Assumptions.** Lifting line, straight-tapered planform, a0 constant along the span, no fuselage; the measured values are those in the shared data file (TN-1422 Table I).

(a)

| a0 | CL_α (per degree), no washout | with 2° washout |
|---|---|---|
| 2π | 0.08889 | 0.08889 |
| 0.95 × 2π | 0.08523 | 0.08523 |
| 0.90 × 2π | 0.08149 | 0.08149 |

(b) The measured 0.085–0.086 per degree corresponds to **a0 ≈ 0.95–0.96 × 2π** (about 6.0 per radian). That is plausible: viscous sections at a few million Reynolds number have lift slopes a few per cent below thin-airfoil theory (a thicker boundary layer on the upper surface de-cambers the section; `P10.06`). Washout does **not** change the slope: twist adds a *basic* load that is independent of α (the problem is linear), so it shifts CL at given α but not dCL/dα.

(c) (1) The section slope a0 at the test Reynolds and Mach numbers is not 2π and is not known exactly. (2) Lifting-line theory ignores the chordwise extent of the tips (lower slope, `P10.11`) and, at AR 9, sweep effects are absent but tip shape matters. (3) The measurement includes tunnel-wall corrections, model support interference and a reading uncertainty; the test model's dihedral (3°) and any fuselage or fairing also change the slope.

**Common mistakes.** Comparing per-radian and per-degree slopes; attributing the whole difference to "experimental error".

## 6. When not to trust the lifting line

(a) At AR 2 the downwash varies strongly along the chord and the span, the tips occupy a large part of the wing, and (for a delta) the wing is highly swept: the 2-D-section-plus-downwash picture fails. The solver will return a finite, smooth but wrong lift slope (too high), with no warning. A vortex lattice (`P10.11`) or slender-wing theory is the right tool.

(b) The lattice: it resolves the chordwise extent of the tips, which lifting-line theory ignores, and its result converges at a measured order. For the trainer's induced drag the difference between e = 0.934 and 0.970 is 4 % of the induced drag, about 4 N at cruise, smaller than the uncertainty of the aircraft's Oswald factor; both are acceptable for performance, the lattice for anything that depends on the tip loading.

(c) The integral equation is smoothing: the load Γ(y) produced by a step in chord or α_L0 is continuous (an infinite step in Γ would need an infinite trailing vortex strength and infinite downwash, which the equation forbids), so only its slope jumps. The Fourier coefficients still converge, though slowly, so a finite answer is found.

## Check-yourself answers

1. The 2-D section lift law at the effective angle, and the downwash of the trailing vortex sheet computed by Biot–Savart.
2. Every sin nθ vanishes at θ = 0 and π, the tips, where Γ must be zero.
3. A1 sets the lift (CL = πAR A1); A2 the rolling moment (Cl = (πAR/4)A2).
4. δ = Σ_{n≥3} n(Aₙ/A1)², e = 1/(1 + δ). δ is a sum of squares with positive weights, so it is never negative: the elliptic load is the minimum.
5. 2π/(1 + 2/6) = 4.712 per radian (0.0822 per degree).
6. The convergence follows the smoothness of the problem in θ: the rectangular chord is smooth across the root, the tapered chord has a slope discontinuity there, which caps the coefficient decay and gives order 2.
7. At the root: its local cl/CL is highest there (1.12 at cruise), because the downwash is smallest at the root.
8. Low aspect ratio (below about 4), swept or delta wings, near or beyond the stall, and wherever the tip or flap-edge detail matters.
