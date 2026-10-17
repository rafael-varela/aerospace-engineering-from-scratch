# P08.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. All results are for inviscid, incompressible, irrotational plane flow.

## 1. A lifting cylinder (s = 15: U = 35 m/s)

**Assumptions.** Potential flow, sea-level density 1.225 kg/m³, Γ clockwise.

(a) Γ = 0.3 × 4π × 35 × 0.1 = **13.195 m²/s**. sin θ_s = −0.3: **θ_s = −17.46° and −162.54°** (both on the lower surface).

(b) Cp = 1 − (2 sin θ + 0.6)². Maximum **1.0** at the two stagnation points; minimum **−5.76 at θ = 90°** (top); at the bottom, θ = 270°, Cp = 1 − (−2 + 0.6)² = **−0.96**. The top is at much lower pressure than the bottom: the clockwise circulation adds to the free-stream speed over the top and subtracts from it underneath, and the net pressure force points up.

(c) Kutta–Joukowski ρUΓ = 1.225 × 35 × 13.195 = **565.72 N/m**; `force_by_pressure` gives **565.72 N/m** of lift and a drag of **10⁻¹⁴ N/m** (round-off). c_l = L′/(½ρU²·2a) = Γ/(Ua) = **3.770 = 4π × 0.3** ✓.

**Independent check.** `force_blasius` on any circle enclosing the cylinder gives the same lift; the stagnation points have zero velocity in `complex_velocity`.

**Common mistakes.** Putting the stagnation points on the upper surface (the sign of Γ); basing c_l on the radius.

## 2. Derivations

(a) On z = ae^{iθ}: a²/z = ae^{−iθ}, so W = U(ae^{iθ} + ae^{−iθ}) = 2Ua cos θ, which is real: **ψ = 0 on the circle**. Far away a²/z → 0 and W → Uz, uniform flow. The vortex term (iΓ/2π) log z = (iΓ/2π)(ln r + iθ) has ψ = (Γ/2π) ln r, constant on r = a: the circle is still a streamline.

(b) dW/dz = U(1 − a²/z²) + iΓ/(2πz). On the circle the velocity is tangential; with u_θ the counterclockwise component, u_r − iu_θ = e^{iθ} dW/dz = U(e^{iθ} − e^{−iθ}) + iΓ/(2πa) = 2iU sin θ + iΓ/(2πa). So u_r = 0 and **u_θ = −2U sin θ − Γ/(2πa)**, and **Cp = 1 − (2 sin θ + Γ/(2πaU))²**. The pressure is p = p∞ + ½ρU²Cp, and F = −∫₀^{2π} p (cos θ, sin θ) a dθ. Expanding Cp, only terms odd in sin θ survive against sin θ, and nothing survives against cos θ (the expression is even about θ = 90°): **D = 0** and L = −½ρU²a ∫ (−4 sin θ · Γ/(2πaU)) sin θ dθ = ½ρU²a × (4Γ/(2πaU)) × π = **ρUΓ**. Stagnation where u_θ = 0: **sin θ_s = −Γ/(4πUa)**; the two points merge at the bottom when **Γ = 4πUa**.

(c) With n ds = −i dz on a counterclockwise contour, F_x + iF_y = −∮ p n ds = i∮ p dz; the constant p∞ + ½ρU² integrates to zero round a closed contour, so F_x + iF_y = −(iρ/2)∮ |dW/dz|² dz. Conjugating, F_x − iF_y = (iρ/2)∮ |dW/dz|² dz̄. On the body, a streamline, dW = (dW/dz) dz is real, so (dW/dz)dz = conj((dW/dz) dz) and |dW/dz|² dz̄ = (dW/dz)² dz. Hence **F_x − iF_y = (iρ/2)∮(dW/dz)² dz**, and since (dW/dz)² is analytic outside the body, the contour can be moved to any enclosing circle. Far away, dW/dz = U + (m + iΓ)/(2πz) + O(z⁻²), with m = 0 for a closed body; (dW/dz)² = U² + iUΓ/(πz) + O(z⁻²). Only the 1/z term contributes: ∮ = 2πi × iUΓ/π = −2UΓ, so F_x − iF_y = (iρ/2)(−2UΓ) = −iρUΓ: **F_x = 0, F_y = ρUΓ** for any shape.

**Common mistakes.** Using the inward normal; forgetting that the conjugate trick needs the contour to be a streamline before it is moved.

## 3. Where to put a static port

(a) ψ = Ur sin θ + mθ/(2π). The stagnation point on the negative x axis has θ = π, y = 0: ψ = m/2. So the body is **ψ = m/2**. Far downstream on the upper surface θ → 0 and Uy → m/2: **b = m/(2U)**. On the axis upstream, u = U + m/(2πx) = 0 at **x = −m/(2πU) = −b/π**.

(b) With U = 1, m = 2π: b = π, D = 2π, stagnation at x = −1.000. The surface is r(θ) = m(π − θ)/(2πU sin θ). The normal velocity on 20 001 traced points, from the finite-difference tangent, is at most **5 × 10⁻⁸ U** (the error of the tangent), so the traced surface carries no flow through it ✓.

(c) Minimum **Cp = −0.587**, at x = 1.04, **0.33 D behind the nose**. |Cp| < 0.01 only beyond x = 200, **32 D from the nose**; |Cp| < 0.001 beyond **307 D**. Far downstream the surface speed is U + m/(2πx), so Cp ≈ −m/(πUx) = −2b/(πx), which reaches −0.01 at x = 200b/π: the half-body's suction dies away only as 1/x.

(d) The half-body is a nose with nothing behind it. A real probe has a stem (or a mast, or the fuselage it is mounted on) downstream, whose own stagnation region raises the pressure ahead of it; the static ports are placed where the nose's suction and the stem's overpressure cancel, a few diameters back. Add the angle of attack, the probe's own boundary layer and the aircraft's local flow field (the static-port position error of `P08.04`), and no simple model gives the reading to the accuracy airworthiness needs: every probe, and every installation, is calibrated.

**Common mistakes.** Taking ψ = 0 as the body (that is the axis upstream); using the principal-branch angle, which puts a cut along the body.

## 4. Your own aeroplane (s = 15: 730 kg, V = 57.5 m/s)

(a) L′ = 730 × 9.80665/10 = **715.89 N/m**; Γ = L′/(ρV) = 715.89/(0.90913 × 57.5) = **13.695 m²/s**; c_l = L′/(½ρV²c) = **0.3176**.

(b) α = asin(0.3176/(2π)) = **2.897°**. The equivalent cylinder (a = 0.375 m) has its stagnation points at **−2.897° and −177.10°**: the rear stagnation point sits exactly at θ = −α. In the cylinder plane with the free stream along x, that is the same as rotating the cylinder flow by α and putting the stagnation point at ζ = a, the Kutta condition. The two descriptions are the same flow.

(c) `force_blasius` on r = 2a: lift **715.885 N/m**, drag 10⁻¹³: equals ρVΓ ✓. In the plate plane (Kutta circulation 13.695 m²/s), body axes: **F_x = −36.18 N/m** along the plate (pointing forward, upstream) and **F_y = 714.97 N/m** normal to it; wind axes: **lift 715.885 N/m, drag 6 × 10⁻¹³ N/m**. The pressure on a plate of zero thickness can only act normal to it, so it would give a force tilted back by α with a drag component; the forward body-axis force is the **leading-edge suction**, the limit of the infinite suction at the sharp edge, which exactly cancels that drag. A real thin aerofoil gets this force from its rounded nose.

**Common mistakes.** Reporting the plate-plane F_x as drag; using c = 2a instead of 4a.

## 5. The Rankine oval family

(a)

| m/(πUa) | l/a | t/a | l/t | max speed/U |
|---|---|---|---|---|
| 0.1 | 1.0488 | 0.1429 | 7.340 | 1.288 |
| 0.5 | 1.2247 | 0.5384 | 2.275 | 1.421 |
| 1 | 1.4142 | 0.8603 | 1.644 | 1.575 |
| 2 | 1.7321 | 1.3065 | 1.326 | 1.739 |
| 5 | 2.4495 | 2.1642 | 1.132 | 1.880 |
| 10 | 3.3166 | 3.1105 | 1.066 | 1.937 |
| 100 | 10.0499 | 9.9834 | 1.007 | 1.993 |

(b) `force_by_pressure` on 256 traced points gives forces between 0 and 2.4 × 10⁻¹⁵ (in ρU²a units) for every oval: **zero**, d'Alembert.

(c) As m grows, the body grows far beyond the source–sink spacing 2a, and seen from its own size the pair is a doublet: the oval becomes a cylinder, whose maximum speed is 2U. As m shrinks the body becomes thin and disturbs the stream little, so its surface speed tends to U (the slender oval at 0.1 still has 1.29U, because its ends are blunt). The slender oval (l/t ≈ 7) is the better fuselage model. Potential flow still predicts zero drag for it: the real drag of a streamlined body is mostly skin friction in the boundary layer (`P08.10`, `P08.14`) plus a little pressure drag from the boundary layer's displacement, all of it viscous.

**Common mistakes.** Expecting the fineness ratio to scale linearly with m; tracing the surface with the principal-branch stream function.

## 6. What potential flow gets right and wrong

(a) Potential flow leaves out **viscosity**: the boundary layer, skin friction and, above all, separation, which leaves a low-pressure wake behind bluff bodies (`P08.16`). It is still useful for wings because, in attached flow, the boundary layer is thin and the pressure outside it is very nearly the potential-flow pressure. Viscosity's main effect on a wing is to fix the circulation through the Kutta condition, which potential flow can impose. Lift and the pressure distribution are predicted well; drag must come from boundary-layer theory.

(b) Γ = 2πa × 2U = 4πUa, so c_l = Γ/(Ua) = **4π = 12.6**, and the two stagnation points just merge at the bottom. The measured lift is far lower because the circulation of a real spinning cylinder is not set by its surface speed. The boundary layer drags a thin layer round, but the flow separates, the wake is shifted, and the effective circulation is a fraction of the "no-slip" value. The failing assumption is that the vortex strength is known from the rotation: in a real fluid it is fixed by viscous separation, which potential flow cannot model.

(c) Trust the **lift** (c_l ≈ 2π sin α, close to measurements in attached flow; thickness and camber shift it) and the **stagnation point** (near the leading edge, slightly below it). Do not trust the **pressure near the leading edge** of a flat plate: the infinite suction is an artefact of zero thickness and a sharp edge; a real nose has a finite suction peak. Do not trust the **drag**: zero in theory, skin friction and form drag in reality.

**Common mistakes.** Concluding that potential flow is "wrong"; it is exact for its assumptions, and the assumptions hold outside thin boundary layers on streamlined bodies.

## Check-yourself answers

1. W(z) = φ + iψ, an analytic function of z = x + iy whose real part is the velocity potential and whose imaginary part is the stream function; dW/dz = u − iv.
2. Laplace's equation (for φ and ψ) and the boundary condition of zero normal flow are linear, so sums of solutions are solutions; the pressure (Bernoulli) is not linear and must be computed from the total velocity.
3. Uniform flow plus a doublet of strength κ = 2πUa² with its axis along the stream (plus a vortex for lift).
4. At sin θ = −Γ/(4πUa) on the surface (two points, both on the lower half for clockwise Γ); at Γ = 4πUa they merge at the bottom; above that a single stagnation point lies in the flow below the cylinder.
5. L′ = ρUΓ per unit span, perpendicular to the free stream, for any closed body in plane potential flow, with Γ the clockwise circulation (aerodynamics convention) for flow from left to right.
6. A closed body in steady inviscid potential flow feels no drag. Potential flow leaves out viscosity: the boundary layer, skin friction and separation, which create real drag.
7. It fixes the circulation so that the flow leaves the sharp trailing edge smoothly (the rear stagnation point at the edge): Γ = 4πUa sin α = πcU sin α, giving c_l = 2π sin α.
8. A body surface is crossed by no flow, and ψ is constant along any curve with no flow across it. The tests trace the surface independently (as points), compute its normal from the traced geometry by FFT differentiation, and check that the velocity from `complex_velocity` has no component along it.
