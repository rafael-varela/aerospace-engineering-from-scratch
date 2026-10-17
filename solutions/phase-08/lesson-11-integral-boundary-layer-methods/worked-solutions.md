# P08.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown.

## 1. The wing's leading edge (s = 15: r = 22.5 mm)

**Assumptions.** Two-dimensional laminar stagnation-point flow; potential-flow edge velocity of a cylinder of the leading-edge radius; ν = 1.8630 × 10⁻⁵ m²/s, μ = ρν = 1.694 × 10⁻⁵ Pa·s.

(a) a = 2U/r = 120/0.0225 = **5333 s⁻¹**. Thwaites at Ue = ax: θ² = 0.45ν∫₀ˣ(ax′)⁵dx′/(ax)⁶ = 0.075ν/a, so **θ = 16.19 μm** (constant along the stagnation region), **λ = 0.075**, and the fits give **H = 2.358, l = 0.3276**, **δ* = 38.17 μm**.

(b) Hiemenz (`falkner_skan(1.0)`): θ√(a/ν) = 0.2923, δ*√(a/ν) = 0.6479, f″(0) = 1.23259, so **θ = 17.28 μm, δ* = 38.29 μm**. Thwaites is **6.3 % low on θ but only 0.3 % low on δ*** because its H(λ = 0.075) = 2.358 is 6 % above the exact 2.216: the two errors cancel in δ* = Hθ. That is luck, not design; the method is built to get θ right on average over many flows and is least accurate in strongly accelerated ones.

(c) τ_w = μUe l/θ = **3.66 Pa** (Thwaites) against μ a x f″(0)√(a/ν) = **3.77 Pa** (exact) at x = 2 mm: Thwaites 2.9 % low. sin(x/r)/(x/r) − 1 = **−0.13 %** at x = 2 mm, so the linear Ue is excellent there.

**Check.** θ = √(0.075ν/a) has units √(m²s⁻¹ · s) = m. A sharper leading edge (smaller r, larger a) gives a thinner layer, θ ∝ √r.

**Common mistakes.** Plugging Ue = 0 into the quadrature and getting 0/0 (or a division error); comparing Thwaites's θ with Hiemenz's θ_η without the √(ν/a) factor.

## 2. Derivations

(a) Integrate u u_x + v u_y = Ue Ue′ + ν u_yy from y = 0 to h > δ. By parts, ∫v u_y dy = Ue v(h) − ∫u v_y dy = Ue v(h) + ∫u u_x dy (continuity, v_y = −u_x; v(0) = 0), and v(h) = −∫₀ʰu_x dy. The viscous term gives ν[u_y]₀ʰ = −τ_w/ρ. So ∫(2u u_x − Ue u_x − Ue Ue′)dy = −τ_w/ρ. Writing 2u u_x − Ue u_x − UeUe′ = −∂/∂x[u(Ue − u)] − Ue′(Ue − u) and integrating:
d/dx(Ue²θ) + Ue Ue′ δ* = τ_w/ρ,  i.e.  **dθ/dx + (2 + H)(θ/Ue) dUe/dx = c_f/2.**
The wall conditions enter through v(0) = 0 and the wall shear; continuity removes v.

(b) Multiply by 2Ueθ/ν: Ue d(θ²)/dx + 2(2 + H)θ²Ue′/ν = c_f Ue θ/ν = 2τ_wθ/(μUe). So **Ue d(θ²)/dx = 2[l − (2 + H)λ] ≡ F(λ)**. With F ≈ 0.45 − 6λ: Ue d(θ²)/dx + 6θ²Ue′ = 0.45ν. Multiply by Ue⁵: d(θ²Ue⁶)/dx = 0.45νUe⁵, and integrating, **θ²Ue⁶ = θ₀²Ue₀⁶ + 0.45ν∫Ue⁵dx**.

(c) Wedge: ∫₀ˣC⁵x′^(5m)dx′ = C⁵x^(5m+1)/(5m + 1), so θ² = 0.45νx/((5m + 1)Ue) and λ = θ²mUe/(xν) = **0.45m/(5m + 1)**, constant. Howarth: ∫₀ˣU₀⁵(1 − x′/L)⁵dx′ = U₀⁵L[1 − (1 − x/L)⁶]/6; with Ue′ = −U₀/L, **λ = −0.075[(1 − x/L)^(−6) − 1]**. Setting λ = −0.09: (1 − x_s/L)^(−6) = 2.2, **x_s/L = 1 − 2.2^(−1/6) = 0.12314**.

(d) With x = Rφ, Ue = 2U sin φ and Ue′ = (2U/R)cos φ: θ² = 0.45νR∫₀^φ(2U)⁵sin⁵φ′dφ′/(2U sin φ)⁶ and λ = θ²Ue′/ν = 0.45 cos φ ∫₀^φ sin⁵φ′dφ′ / sin⁶φ. U, R and ν cancel: λ depends on φ only, so the separation angle (**103.1°** by the code) is universal for a laminar layer in the potential-flow pressure field. The test `test_cylinder_separation_is_independent_of_speed_size_and_viscosity` checks it to 10⁻⁶.

**Common mistakes.** Dropping the UeUe′δ* term (it is what makes H appear); taking 0.45 − 6λ as exact (it is a fit to exact solutions).

## 3. A cabin air-duct diffuser (s = 15: l₀ = 0.35 m)

**Assumptions.** Laminar, two-dimensional wall layer (the duct is wide compared with δ), the core velocity imposed as given, sea-level ν.

(a) ∫₀^ξU₀⁵(1 + ξ′/L)^(−5)dξ′ = U₀⁵L[1 − (1 + ξ/L)^(−4)]/4; θ² = 0.1125νL(1 + ξ/L)⁶[1 − (1 + ξ/L)^(−4)]/U₀ and Ue′ = −U₀/(L(1 + ξ/L)²), so **λ = −0.1125[(1 + ξ/L)⁴ − 1]**. λ = −0.09 at (1 + ξ/L)⁴ = 1.8: **ξ_s = 0.1583 L = 0.475 m**, where **Ue has fallen 13.7 %**.

(b) θ₀ = √(0.45νl₀/U₀) = **0.480 mm**; λ at the entrance = θ₀²Ue′/ν = −0.45 l₀/L_d = **−0.0525**. Separation at **ξ_s = 0.158 m**, after Ue has fallen only **5.0 %**. A layer that arrives already thick starts the diffuser more than halfway to the separation value of λ, because λ ∝ θ².

(c) Marching the composite Ue from 10⁻⁴ m: separation **0.153 m** after the diffuser entrance (Ue down 4.8 %), against Thwaites's 0.158 m over the same composite: **3 % apart**, the usual gap between the closure and the PDE.

(d) Shorten or remove the straight section ahead of the diffuser, reduce its angle (larger L_d), or make the layer turbulent before the diffuser (a trip), since a turbulent layer tolerates far more deceleration. A separated layer blocks part of the cross-section, so the core does not slow as designed and the static-pressure recovery falls; the flow may also become unsteady and noisy.

**Common mistakes.** Using the Blasius θ₀ instead of Thwaites's (it changes ξ_s by 3 %, 0.163 m against 0.158 m; either is defensible if stated); applying Thwaites after λ has passed −0.09.

## 4. How good is the closure? (s = 15: m₁ = 0.07)

| m | β | θ√Re_x/x exact / Thwaites | c_f√Re_x exact / Thwaites | H exact / Thwaites |
|---|---|---|---|---|
| 0.07 | 0.1308 | 0.5830 / 0.5774 (−0.97 %) | 0.9060 / 0.8856 (−2.25 %) | 2.456 / 2.525 |
| −0.05 | −0.1053 | 0.7515 / 0.7746 (+3.08 %) | 0.4270 / 0.4413 (+3.36 %) | 2.818 / 2.753 |
| 0.5 | 0.6667 | 0.3779 / 0.3586 (−5.11 %) | 1.7994 / 1.7486 (−2.83 %) | 2.262 / 2.391 |

(b) The errors are smallest near the flat plate (θ +1.0 %, c_f −1.2 % at m = 0) and grow towards both ends: about −5 % in θ for strongly accelerated flow (m = 0.5; −6.3 % at Hiemenz, m = 1) and +3 to +4 % in θ for decelerating flow approaching separation (+4.4 % at m = −0.0904). A one-parameter closure cannot be exact everywhere because real profiles at the same λ are not all the same shape.

(c) The value −0.068 belongs to a Falkner–Skan profile that has been in the same equilibrium for its whole history. Howarth's layer reaches separation with a different history (it started as a flat-plate layer and was decelerated progressively), so its profile at a given λ is not the wedge profile at that λ. The value −0.09 is an empirical criterion that works well for such progressively decelerated layers (2.7 % late for Howarth's flow, 1.4° early on the cylinder against the marched PDE), even though it is not the λ at which a wedge profile loses its shear. An empirical criterion is right for the kind of flow it was found to work on: trust it in flows like those, and check it by an independent route elsewhere.

**Common mistakes.** Comparing θ_η of Falkner–Skan directly with Blasius's 0.664 (different scalings: √2 apart at m = 0).

## 5. Two routes to separation

(a) x_s/L = **0.12047, 0.12036, 0.12024, 0.12009** for nx = 1000–8000. The values fall slowly and steadily (the x-error is first order but the extrapolation to τ = 0 also depends on the last few stations), heading for **about 0.1199–0.1200**, the published 0.1199. Thwaites gives **0.1231**, 2.7 % late.

(b) **104.84°, 104.69°, 104.61°, 104.54°** for the four grids, converging towards about 104.5°; Thwaites gives **103.1°**, about 1.4° early. (Both are for the *potential-flow* pressure distribution; the real cylinder separates far earlier, because its wake changes the pressure distribution: the outer flow and the layer interact.)

(c) The boundary-layer equations are parabolic only while u > 0: information is carried downstream by the flow. With reversed flow near the wall, information travels upstream too, so a downstream march is ill-posed and the implicit solve fails (Goldstein's singularity: τ_w ∝ √(x_s − x) with infinite dτ_w/dx at x_s). Going further needs either the full Navier–Stokes equations (`P08.18`) or an inverse boundary-layer method that prescribes δ* and computes Ue, coupled to the outer flow (viscous–inviscid interaction, `P10.06`).

**Common mistakes.** Reading the step at which marching stopped as the separation point (it is a grid-dependent overshoot or undershoot; extrapolate τ² to zero).

## 6. A laminar-flow wing

(a) Just behind the peak λ = θ_peak² (dUe/dx)/ν. Moving the peak from 0.5c to 0.7c makes the layer at the peak thicker (θ = 0.236 → 0.280 mm, so θ² rises 40 %) *and* squeezes the same 0.3U of deceleration into 0.45 m instead of 0.75 m (dUe/dx from −24 to −40 s⁻¹). λ jumps from −0.072 to −0.168, past −0.09 at once. A laminar-flow aerofoil buys a long run of laminar flow at the price of a steep pressure recovery that only a turbulent layer can survive.

(b) The separated laminar shear layer is very unstable; it becomes turbulent within a short distance and the turbulent layer usually reattaches, enclosing a **laminar separation bubble**. Bubbles cost pressure drag and some friction, and at low Reynolds number (gliders, drones; `P10.14`) a long bubble can burst and stall the section. Designers place a trip or shape the recovery so that transition happens just ahead of laminar separation.

(c) The boundary-layer approximation itself: near separation the layer thickens rapidly, ∂p/∂y is no longer negligible and the outer flow is changed by the layer (the pressure is no longer imposed). On a model, a row of surface pressure taps shows it: a laminar separation bubble appears as a flat pressure plateau followed by a sharp recovery at reattachment, not the smooth distribution assumed.

## Check-yourself answers

1. λ = θ²(dUe/dx)/ν, the ratio of the pressure-gradient force to the viscous force across the layer. λ > 0: accelerating (favourable) flow, thin full profiles; λ < 0: decelerating (adverse) flow; λ = −0.09: Thwaites's separation.
2. dθ/dx + (2 + H)(θ/Ue)dUe/dx = c_f/2: growth of the momentum deficit, plus the effect of the pressure gradient on the deficit (2θ from momentum, Hθ = δ* from the pressure force), equals the wall shear.
3. That l and H are universal functions of λ alone (one parameter describes every profile), and that F(λ) = 2[l − (2 + H)λ] is the straight line 0.45 − 6λ.
4. θ is 1.0 % high (0.6708 against 0.6641 x/√Re_x), c_f 1.2 % low (0.6559 against 0.6641/√Re_x) and H = 2.61 against 2.591.
5. Ue = ax grows as fast as diffusion thickens the layer: the inflow towards the wall balances viscous growth, so θ² = 0.075ν/a is constant.
6. From a fresh start, about 12 % (Howarth: published 12.0 %, Thwaites 12.3 %; the diffuser of Problem 3, 13.7 %). A layer that arrives with a finite θ starts closer to λ = −0.09, because λ ∝ θ²: 5 % in Problem 3.
7. With reversed flow near the wall the equations are no longer parabolic in x; information travels upstream, and a downstream march has nothing to march with (Goldstein singularity at separation).
8. Thwaites and the marched PDE give 103.1° and about 104.5° for the potential-flow pressure distribution. A real subcritical cylinder separates much earlier because its separated wake changes the surface pressure from the potential-flow one; the layer and the outer flow must be solved together.
