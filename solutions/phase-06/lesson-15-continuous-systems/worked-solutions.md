# P06.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. All results are outputs of the model, not measurements. Material values are typical, not design allowables.

## 1. An instrument mast (s = 15: L = 2.5 m)

**Assumptions.** Uniform Euler–Bernoulli tube, rigid clamp, point-mass sensor without rotary inertia.

(a) A = π(0.06² − 0.056²)/4 = **3.644 × 10⁻⁴ m²**; I = π(0.06⁴ − 0.056⁴)/64 = **1.534 × 10⁻⁷ m⁴**; EI = **15 342 N·m²**; ρA = **0.5831 kg/m**; mast mass **1.458 kg**; μ = 0.2/1.458 = **0.1372**.

(b) Bare mast: **f = 14.52, 91.02, 254.85 Hz**.

(c) With the sensor: **f = 11.65, 77.80, 225.47 Hz**.

(d) SDOF estimate **11.716 Hz**, **+0.59 %**: between the lesson's values at μ = 0.1 (+0.73 %) and μ = 0.25 (+0.35 %), as it should be.

**Check.** f₂/f₁ = 6.27 for the bare mast, the cantilever's (4.694/1.875)².

**Common mistakes.** Using D/2 in I as if it were a diameter formula; forgetting that the mass ratio, not the sensor mass alone, sets the effect.

## 2. The beam toolkit

(a) On a slice dx with shear V and moment M (sign conventions of `P05.12`), vertical equilibrium with d'Alembert's force: dV/dx = q − ρAẅ; moments: dM/dx = V; and M = EIw″. So (EIw″)″ = q − ρAẅ, which with q = 0 is **(EIw″)″ + ρAẅ = 0**.

(b) φ = C₁ cosh βx + C₂ sinh βx + C₃ cos βx + C₄ sin βx. φ(0) = 0: C₁ + C₃ = 0. φ′(0) = 0: β(C₂ + C₄) = 0. So φ = C₁(cosh − cos) + C₂(sinh − sin). φ″(L) = 0: C₁(cosh + cos) + C₂(sinh + sin) = 0; φ‴(L) = 0: C₁(sinh − sin) + C₂(cosh + cos) = 0 (all at βL). Determinant: (cosh + cos)² − (sinh² − sin²) = cosh² − sinh² + cos² + sin² + 2 cos cosh = 2(1 + cos βL cosh βL) = 0: **(4)**. From the first row, C₂/C₁ = −(cosh + cos)/(sinh + sin) = −σ: the mode shape.

(c) The tip mass m_t moves with the beam's end, w(L, t). The only force on it is the beam's shear at the end, so m_t ẅ(L) = EI w‴(L) with the free-end sign convention (positive shear pushing the mass up), i.e. **EI φ‴(L) = −m_t ω² φ(L)** with w = φ cos ωt. The moment condition is unchanged. In dimensionless form the fourth row becomes φ‴(L) + μβ⁴L φ(L) = 0 (using ω² = β⁴EI/ρA), adding μβL(cos sinh − sin cosh) terms to the determinant: **(5)**. μ = 0 recovers (4).

(d) ∫₀ᴸ φⱼφᵢ⁗ dx = [φⱼφᵢ‴]₀ᴸ − [φⱼ′φᵢ″]₀ᴸ + ∫φⱼ″φᵢ″ dx. At x = 0 the brackets vanish because φⱼ(0) = φⱼ′(0) = 0 (clamp); at x = L because φᵢ‴(L) = φᵢ″(L) = 0 (free end). So βᵢ⁴∫φⱼφᵢ = ∫φⱼ″φᵢ″, symmetric in i, j; subtracting the swapped equation gives **(βᵢ⁴ − βⱼ⁴)∫φᵢφⱼ dx = 0**. With a tip mass the boundary term at L no longer vanishes, and the orthogonality weight acquires m_t φᵢ(L)φⱼ(L).

(e) With θ = jπh/L and wᵢ = sin(iθ), the stencil gives sin((i−2)θ) − 4 sin((i−1)θ) + 6 sin(iθ) − 4 sin((i+1)θ) + sin((i+2)θ) = sin(iθ)(2 cos 2θ − 8 cos θ + 6) = sin(iθ)(2 − 2 cos θ)². At the ends, w₀ = sin 0 = 0 and the ghost value w₋₁ = −w₁ is exactly sin(−θ): the sine satisfies the boundary rows automatically. Hence **λⱼ = (2 − 2 cos(jπh/L))²/h⁴**, (7).

**Common mistakes.** Sign of the tip-mass shear; forgetting that the ghost-node rule is what makes the sine exact.

## 3. A rocket in flight (s = 15: L = 37.5 m)

(a) Free–free: φ″ = φ‴ = 0 at both ends. The 4 × 4 determinant reduces to **1 − cos βL cosh βL = 0**, the clamped–clamped equation (the boundary conditions are "dual"). The free–free beam also has β = 0 solutions: **rigid translation and rigid rotation** (φ = 1 and φ = x), with zero frequency.

(b) √(EI/(ρAL⁴)) with βL = 4.730, 7.853, 10.996: **f = 2.002, 5.518, 10.818 Hz**.

(c) 60 elements, no supports: the first two eigenvalues are **−2.1 × 10⁻⁷ and 1.3 × 10⁻⁷ s⁻²** (zero within round-off: the two rigid-body modes), then **2.00183, 5.51812, 10.8177 Hz**, matching (b) to 10⁻⁶.

(d) The first mode's deflection is zero at **x/L = 0.224 and 0.776** (closed form 0.22416, 0.77584; FEM 0.22423, 0.77577). An *accelerometer* there sees no mode-1 translation. A *rate gyro* senses rotation, w′, which is largest near the deflection nodes and zero at the deflection antinodes; for mode 1 that is the middle of the vehicle (x/L = 0.5, by symmetry). Real vehicles place gyros where the slope of the troublesome modes is small, or filter the bending frequencies out of the control loop.

**Common mistakes.** Supporting the FE model to avoid a singular K (it then is no longer free–free); placing a rate sensor at a deflection node.

## 4. A model fit for purpose (s = 15: tank 35 kg)

(a) Relative frequency errors (FE − exact)/exact:

| elements | mode 1 | mode 2 | mode 3 | mode 4 | mode 5 |
|---|---|---|---|---|---|
| 3 | 1.0 × 10⁻⁴ | 3.3 × 10⁻³ | 1.2 × 10⁻² | 0.16 | 0.32 |
| 4 | 3.3 × 10⁻⁵ | 1.2 × 10⁻³ | 7.7 × 10⁻³ | 1.5 × 10⁻² | 0.14 |
| 6 | 6.6 × 10⁻⁶ | 2.5 × 10⁻⁴ | 1.8 × 10⁻³ | 6.4 × 10⁻³ | 1.5 × 10⁻² |
| 8 | 2.1 × 10⁻⁶ | 8.0 × 10⁻⁵ | 6.1 × 10⁻⁴ | 2.2 × 10⁻³ | 5.8 × 10⁻³ |
| 10 | 9 × 10⁻⁷ | 3.3 × 10⁻⁵ | 2.5 × 10⁻⁴ | 9.5 × 10⁻⁴ | 2.5 × 10⁻³ |
| 12 | 4 × 10⁻⁷ | 1.6 × 10⁻⁵ | 1.2 × 10⁻⁴ | 4.7 × 10⁻⁴ | 1.3 × 10⁻³ |
| 16 | 1 × 10⁻⁷ | 5.1 × 10⁻⁶ | 4.0 × 10⁻⁵ | 1.5 × 10⁻⁴ | 4.1 × 10⁻⁴ |
| 20 | 1 × 10⁻⁷ | 2.1 × 10⁻⁶ | 1.6 × 10⁻⁵ | 6.3 × 10⁻⁵ | 1.7 × 10⁻⁴ |

10⁻⁴ needs about 4, 8, 16 and 20 elements for modes 1–4 (3 elements leave mode 1 at 1.01 × 10⁻⁴); mode 5 is still at 1.7 × 10⁻⁴ with 20. Rule: **four to five elements per mode number** for 10⁻⁴ (roughly two to three elements per half-wavelength at the highest mode wanted), and a model has no useful information about modes beyond about its number of elements divided by two.

(b) Baseline 6.224 and 39.004 Hz. The 35 kg tank lowers f₂ most at **x = 2.1 m** (f₂ = 25.58 Hz, f₁ = 5.68 Hz). In the outboard half f₂ is left essentially unchanged at **x = 3.9 m** (f₂ = 38.998 Hz), where f₁ falls to 3.92 Hz. The node of mode 2 from `mode_shape` is at **x = 3.917 m**: the same place, within the 0.05 m grid.

(c) Treat the mass as a small change: the Rayleigh quotient of mode i becomes ωᵢ² ≈ (stiffness term)/(1 + m φᵢ(x)²/Mᵢ), with Mᵢ the modal mass. The drop in ωᵢ² is proportional to φᵢ(x)²: largest at the mode's largest inboard deflection (mode 2 near 2 m) and zero at its node (3.92 m). Mode 1 has no node and its deflection grows towards the tip, so f₁ falls more the further out the tank goes.

**Common mistakes.** Using the same element count for every mode; looking for the mode-2 node with a lumped model too coarse to place it.

## 5. Is Euler–Bernoulli enough?

(a) Mode n of a cantilever has a half-wavelength of roughly L/(n − ½) (from βL ≈ (n − ½)π). Ten depths is 1.8 m, so n − ½ ≈ 2.8: about **mode 3**. Beyond that, the shear strain and the rotational inertia of the cross-sections, both ignored by Euler–Bernoulli, become a few per cent of the bending terms, and the theory overestimates the frequencies (Timoshenko theory, Extension 1).

(b) (1) **The skins and stringers** add bending stiffness: raises the frequency (f ∝ √EI; 8.5/6.22 needs EI × 1.87). (2) **The mass distribution**: if the real wing's mass sits further inboard than the uniform 9 kg/m, f₁ rises (experiment 4). (3) **The root**: a flexible root would *lower* f₁, so it cannot explain a higher measurement; and coupling with torsion shifts the measured mode too. The model's spar-only stiffness is the first suspect.

(c) From experiment 1, 8 lumped elements give f₁ **0.7 % low** and f₃ **3.9 % low**, both from below (consistent mass at 8 elements: +2 × 10⁻⁶ and +6 × 10⁻⁴, from above). Acceptable for a first estimate of f₁, not for f₃.

(d) Trust: the scaling laws (f ∝ √(EI/ρA)/L²), the relative effects (where a mass moves which mode, the nodes), and the numerical accuracy of the FE model, all of which this lesson verified independently. Do not trust as design values the absolute frequencies of the real wing: the model's EI is the spar alone, and bending–torsion coupling and root fixity are missing. Those need a ground vibration test.

## Check-yourself answers

1. EI w⁗ + ρA ẅ = 0; both terms in N/m (force per unit length).
2. Clamped: w = 0, w′ = 0. Pinned: w = 0, w″ = 0. Free: w″ = 0, w‴ = 0. A tip mass changes the shear condition to EI w‴(L) = m_t ẅ(L).
3. 1 + cos βL cosh βL = 0. As βL grows, cosh βL grows exponentially, so cos βL must approach zero: βL → (2n − 1)π/2.
4. ω ∝ √(EI/ρA)/L²: as L⁻², √EI and 1/√(ρA).
5. The operator d⁴/dx⁴ with the standard boundary conditions is self-adjoint (integration by parts twice leaves symmetric boundary terms that vanish), so modes with different βᵢ are orthogonal with weight ρA (plus any point masses).
6. The bending energy contains w″, so the trial functions must be continuous with continuous slopes (C¹); a Hermite element uses the deflection and the slope at each of its two nodes.
7. Consistent mass keeps the FE method a Rayleigh–Ritz minimisation over admissible functions, so frequencies come from above, at order 4 for cubics. Lumped mass underestimates the kinetic energy's coupling and removes the rotary terms, giving frequencies from below at order 2.
8. Near the second mode's largest inboard deflection (about 0.4 L for a cantilever, 2 m on the 5 m wing), where the first mode's deflection is still small (0.23 of its tip value): a mass there lowers f₂ by about a third while lowering f₁ by under 10 %.
