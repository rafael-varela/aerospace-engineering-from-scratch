# P05.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Material values are typical, not design allowables.

## 1. A bracket arm on edge and on its side (s = 15: h = 45 mm)

**Assumptions.** Euler–Bernoulli cantilever, linear elastic, the 800 N load is the limit load, the support is rigid.

(a) FBD: the arm with the tip load P = 800 N downward and, at the support, an upward reaction of 800 N and a reaction moment. Moment equilibrium about the support: **M = PL = 800 × 0.25 = 200 N·m** (hogging: the top is stretched).

(b) I = bh³/12 = 20 × 45³/12 = **151 875 mm⁴**; S = I/(h/2) = bh²/6 = **6750 mm³**; σ = M/S = 200/6.75 × 10⁻⁶ = **29.6 MPa**, tension on the **top** edge, compression on the bottom.

(c) MS_yield = 345/29.6 − 1 = **+10.6**; MS_ultimate = 483/(1.5 × 29.6) − 1 = **+9.87**. Very large margins: such an arm is usually sized by stiffness or by the fastener holes, not by bending strength.

(d) On its side the depth is 20 mm and the width 45 mm: I = 45 × 20³/12 = 30 000 mm⁴, σ = 200 × 0.010/3 × 10⁻⁸ = **66.7 MPa**, a factor of h/b = **2.25** higher; MS_yield = **+4.18**.

**Independent check.** Units: N·m/m³ = Pa. The ratio of section moduli of a rectangle in its two orientations is h/b, which the numbers reproduce (66.7/29.6 = 2.25).

**Common mistakes.** Using I instead of S (forgetting c); taking c = h instead of h/2; swapping b and h in bh³/12 when the bar is turned.

## 2. Derivations

(a) Compatibility: plane sections give ε = −κy (geometry only). Constitutive law: σ = Eε = −Eκy. Equilibrium, axial: ∫σ dA = −Eκ∫y dA = 0 for pure bending, so **∫y dA = 0: y is measured from the centroid**. Equilibrium, moment: M = −∫σy dA = Eκ∫y² dA = EκI, so **κ = M/(EI)** and σ = −My/I.

(b) With σ = E(a + bx + cy) and x, y from the centroid (∫x dA = ∫y dA = 0):
N = EaA = 0 ⇒ a = 0;
M = −∫σy dA = −E(bI_xy + cI_xx);
M_h = −∫σx dA = −E(bI_yy + cI_xy).
Solving the 2 × 2 system with D = I_xxI_yy − I_xy²:
Ec = −(M I_yy − M_h I_xy)/D, Eb = −(M_h I_xx − M I_xy)/D, so
**σ = −[(M I_yy − M_h I_xy) y + (M_h I_xx − M I_xy) x]/D.**
With I_xy = 0, M_h = 0: σ = −M I_yy y/(I_xx I_yy) = −My/I_xx. ✓ The neutral axis σ = 0 is the line (M I_yy − M_h I_xy) y = −(M_h I_xx − M I_xy) x, at angle **tan α = −(M_h I_xx − M I_xy)/(M I_yy − M_h I_xy)**; for M_h = 0, tan α = I_xy/I_yy.

(c) Each cap: own I = bt³/12 = At²/12, parallel-axis term A(h_c/2)². I = A h_c²/2 + A t²/6. σ_outer = M(h/2)/I. Ratio:

(M/(A h_c)) / (M h/(2I)) = 2I/(A h_c h) = **(h_c² + t²/3)/(h_c h)**, h = h_c + t.

For h_c = 170 mm: t = 10 mm **0.9455**, 5 mm **0.9725**, 2 mm **0.9889** (the `section_properties` values agree to four figures). As t → 0, h → h_c and the ratio → **1**: the idealisation is exact for thin caps with no web. Without a web it is *unconservative* by 1 − ratio (the outer fibre is farther out than the cap centroid); with the lesson's 3 mm web it became conservative because the web's share of I outweighs that.

**Common mistakes.** Forgetting that a ≠ 0 if the origin is not the centroid; sign errors from mixing the stress-resultant convention (M = −∫σy dA here) with a textbook's; dropping the cap's own t²/12 term and then finding a ratio of h_c/h.

## 3. A launch-vehicle tank at maximum dynamic pressure (s = 15: t = 5.5 mm)

**Assumptions.** Thin wall (t/R = 0.003), linear elastic, plane sections, no pressure, the given loads act together.

(a) With y = R sin φ, dA = Rt dφ: A = 2πRt, I = ∫₀^{2π} R² sin² φ · Rt dφ = **πR³t**. A = **622.0 cm²**, I = **0.10077 m⁴**; the exact annulus gives the same to **−2.3 × 10⁻⁶** relative (the error is of order (t/R)²).

(b) σ_N = N/A = −1.5 × 10⁶/0.06220 = **−24.1 MPa** (everywhere). σ_b = MR/I = 2.0 × 10⁶ × 1.8/0.10077 = **±35.7 MPa** at the two points in the plane of bending. Totals: **−59.8 MPa** (compression, on the side the bending compresses, which is the downwind side of the bend) and **+11.6 MPa** (tension, diametrically opposite).

(c) No tension when MR/I ≤ |N|/A, i.e. M ≤ |N| I/(AR) = |N| (πR³t)/(2πRt · R) = **|N|R/2 = 1.350 MN·m**, independent of t (both A and I are proportional to t).

(d) A thin cylindrical shell buckles in compression at a small fraction of the material's yield stress, so the compressive stress is limited by stability, not strength.

**Independent check.** Integrating the thin-wall stress σ(φ) = N/A − (M R sin φ)/I times y dA around the circumference returns M, and times dA returns N.

**Common mistakes.** Using the solid-circle I = πR⁴/4; adding the axial and bending stresses with the wrong sign on one side; treating the tension side as critical.

## 4. Size the spar of your own aeroplane (s = 15: m = 730 kg, semi-span 5.25 m, h = 0.170 m)

**Assumptions.** Elliptic lift, n = 3.8, tail load and wing weight ignored, the spar alone carries bending.

(a) L = 3.8 × 730 × 9.80665/2 = **13 601.8 N**; M₀ = L · 4s/(3π) = **30 307.2 N·m**.

(b) `brentq` on 1.5 M₀ (h/2)/I(b_f) − 572 × 10⁶ gives **b_f = 46.1 mm**, A_cap = **461 mm²**; the two-flange estimate is **497 mm²** (8 % more, conservative).

(c) At MS_ultimate = 0 the limit stress is F_tu/1.5 = **381.3 MPa** by construction, whatever the load; MS_yield = 503/381.3 − 1 = **+0.32**. The seed changes the load and the depth, and the solver changes the cap until the stress is F_tu/1.5 again.

(d) With 20 elements, the first element midpoint stress is **359.36 MPa** from E κ c and **359.23 MPa** from the flexure formula, relative difference **3.6 × 10⁻⁴**, consistent with the lesson's convergence table.

**Independent check.** The fibre model at the solved width returns the same outer-fibre stress to 10⁻⁴.

**Common mistakes.** A `brentq` bracket that does not change sign; comparing the finite-element midpoint stress with the root stress instead of the stress at the same station; forgetting that the solved section's limit stress is fixed.

## 5. A Z-stringer, free and riveted (s = 15: M = 55 N·m)

(a) I_xx = **12 255 mm⁴**, I_yy = **3064 mm⁴**, I_xy = **4567 mm⁴**; free neutral axis at **56.1°** to the horizontal.

(b) Stresses (MPa):

| point (x, y) mm | free (general formula) | symmetric −My/I_xx |
|---|---|---|
| top-flange tip (14.2, 14.2) | **+70.3** | **−63.7** |
| bottom-flange tip (−14.2, −14.2) | **−70.3** | **+63.7** |
| top web corner (−0.8, 14.2) | **−155.4** | **−63.7** |
| bottom web corner (0.8, −14.2) | **+155.4** | **+63.7** |

The free stringer's largest stress is at the web corners and is 2.4 times the symmetric estimate; the flange tips change sign.

(c) With ε = cy: M = −∫Ecy · y dA = −EcI_xx ⇒ σ = Ecy = **−My/I_xx** exactly. The lateral resultant is M_h = −∫Ecy · x dA = −EcI_xy = **M I_xy/I_xx = 20.49 N·m**, which the skin and rivets must supply. `bending_stress_general(55, 20.49, 0.0142, 0.0142, ...)` returns −63.7 MPa, the symmetric value. ✓

(d) A thin skin between widely spaced ribs restrains the stringer only partly; the truth lies between the two analyses. Use the free (general) analysis as the conservative bound for stress, check that the rivets can carry the lateral load implied by M_h, and if the difference matters, model the skin's lateral stiffness (or test a panel).

**Common mistakes.** Using the symmetric formula for a free Z; forgetting that restraint needs a force (M_h) that something must carry.

## 6. Is the spar margin real?

Baseline: elliptic lift, M₀ = 29 655 N·m, σ = 275.0 MPa, MS_ultimate = +0.386.

| assumption removed | effect on root cap stress | number | MS_ultimate |
|---|---|---|---|
| (a) uniform instead of elliptic lift | raises | M = Ls/2 = 34 936 N·m, σ = 324.0 MPa | **+0.18** |
| (b) wing weight included (45 kg at 3.8 g) | lowers | relief qs²/2 = 4192 N·m, σ = 236.2 MPa | **+0.61** |
| (c) 10 % tail download | raises | M = 32 620 N·m, σ = 302.5 MPa | **+0.26** |
| (d) two 6 mm holes in the bottom cap | raises (tension side) | centroid moves up 6.54 mm, I = 8.769 × 10⁶ mm⁴, σ = 326.5 MPa | **+0.17** |
| (e) skins share the bending | lowers | not computed: needs the wing-box section | — |
| (f) compression cap can buckle | may govern | not computable here: `P05.15` and plate buckling | unknown |

Ranking for the decision: (f) first, because it is the only one that could make the margin negative and it is not even estimated; then (d) and (a), each worth about 0.2 of margin; then (c). (b) and (e) help but should not be counted on until computed with the real wing mass and the real wing box. Recommendation: **do not remove material.** Ask first for a cap crippling and web buckling check and a net-section check at the fastener holes, and use a lifting-line (or measured) load distribution with the tail load.

**Common mistakes.** Treating every effect as a reduction; adding margins instead of recomputing stress; ignoring that holes in the tension cap move the neutral axis.

## Check-yourself answers

1. From the centroid of the section (the neutral axis passes through it), fixed by ∫σ dA = 0: pure bending has no net axial force.
2. S = I/c = 9.7 × 10⁶/90 = 1.08 × 10⁵ mm³; M = σS = 100 × 10⁶ × 1.08 × 10⁻⁴ = 10.8 kN·m.
3. I weights area by the square of its distance from the neutral axis; the I-section puts most of its area in the caps, far from the axis, where the rectangle wastes area near it.
4. The top cap is in compression; the moment is sagging (positive) in this lesson's convention.
5. At 4s/(3π) = 0.424 s from the root, against 0.5 s for a uniform load: the elliptic root moment is 15 % smaller for the same lift (4/(3π)/0.5 = 0.849).
6. MS = allowable/(factor × applied) − 1, required ≥ 0; factor 1.0 against yield at limit load and 1.5 against ultimate strength (ultimate load = 1.5 × limit).
7. When the section has a nonzero product of inertia I_xy about its centroidal axes (it is not symmetric about either axis), unless something restrains lateral bending.
8. The two-flange estimate (290.7 MPa) for a first sizing: it is simple, slightly conservative for a spar with a real web, and its error is known; then check with the full flexure formula.
