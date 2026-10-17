# P05.08  Worked solutions

Each solution states assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. A tail-rotor drive tube (s = 15: 175 kW)

**Assumptions.** Steady torque; circular tube; typical 2024-T3 (G = 27.5 GPa, Fty = 345 MPa).

(a) T = 175 000/(6000 × 2π/60) = **278.52 N·m**. J = π(50⁴ − 47⁴)/32 = 1.3453 × 10⁵ mm⁴. τ_outer = 278 520 × 25/J = **51.76 MPa**; τ_inner = 278 520 × 23.5/J = **48.65 MPa**. Twist over 2 m: **8.633°**.

(b) Fsy = 345/√3 = 199.2 MPa; ratio **3.85**.

(c) Bredt–Batho with r_m = 24.25 mm: τ = T/(2πr_m²t) = **50.25 MPa**, between the inner and outer values (it is their average across the wall), 2.9 % below the outer surface stress. Report the outer-surface value from (1) for strength (the material yields where the stress is highest); J itself differs by only 0.1 % (J_Bredt/J = 0.99904).

**Check.** τ_outer/τ_inner = 25/23.5, the radius ratio, as (1) requires.

**Common mistakes.** Forgetting to convert rpm to rad/s; using the mid-line stress for a strength check on a tube that is not very thin.

## 2. Circular shafts, thin closed sections and open sections

(a) **Compatibility:** a circular section rotates rigidly, so a surface line at radius r tilts by γ = r dφ/dx. **Hooke's law:** τ = Gγ = Gr dφ/dx. **Equilibrium:** T = ∫τr dA = G(dφ/dx)∫r² dA = GJ dφ/dx. Hence dφ/dx = T/(GJ) and **τ = Tr/J**, and for a uniform shaft **φ = TL/(GJ)**.

(b) An element of wall between perimeter stations 1 and 2 and length dx along the axis: the only axial forces are the complementary shear flows on the two perimeter cuts, q₁dx and q₂dx, so **q₁ = q₂ = q**. The wall element ds carries force q ds along the wall; about any point O its moment is q ds p = q × 2 dA_triangle; around the cell the triangles sum to A_m, so **T = 2A_m q**. Strain energy per unit length of tube: ∮(τ²/2G) t ds = (q²/2G)∮ds/t. The work of the torque per unit length as the twist rises from 0 is ½T dθ/dx. Equating: dθ/dx = (q²/G)∮ds/t / T = **T∮(ds/t)/(4A_m²G)**, and J = T/(G dθ/dx) = 4A_m²/∮ds/t.

(c) r_o⁴ − r_i⁴ with r_{o,i} = r_m ± t/2: the even powers of t/2 cancel and (r_m + a)⁴ − (r_m − a)⁴ = 8r_m³a + 8r_m a³ with a = t/2, i.e. 4r_m³t + r_m t³. So J_exact = (π/2)(4r_m³t + r_m t³) = 2πr_m³t(1 + t²/(4r_m²)) and **J_Bredt/J_exact = 1/(1 + (t/2r_m)²)**, exactly.

(d) Slit tube: one strip of width 2πr_m and thickness t, J_open = 2πr_m t³/3. Closed: J = 2πr_m³t. Ratio **3(r_m/t)²**.

**Check.** The tests confirm (b) by taking moments about an arbitrary point, (c) to 10⁻¹², and (d) exactly.

**Common mistakes.** Taking moments of the shear stress instead of the shear flow (forgetting the thickness); assuming q varies with thickness; using the outer area instead of the mid-line area A_m.

## 3. The D-spar of a rotor blade (s = 15: T = 550 N·m)

**Assumptions.** Thin walls; free warping; typical 7075-T6 (G = 26.96 GPa).

(a) With a 400-segment nose: **A_m = 2513.2 mm²**, **q = 109.42 N/mm**, **τ_nose = 91.2 MPa**, **τ_web = 54.7 MPa**, **J = 1.7458 × 10⁵ mm⁴**, twist **6.696°/m**.

(b) Exact: A_m = π40²/2 = 2513.3 mm²; ∮ds/t = π40/1.2 + 80/2.0; J = 1.7459 × 10⁵ mm⁴. The polygon is within 10⁻⁴.

(c) Open: J = (π40 × 1.2³ + 80 × 2.0³)/3 = **285.7 mm⁴**; the closed D-spar is **611 times** stiffer.

**Check.** q × 2A_m = 109.42 × 5026.5 = 550 000 N·mm.

**Common mistakes.** Using the polygon's chord length for the web but giving it the nose thickness; forgetting that the shear flow, not the stress, is equal in the nose and the web.

## 4. Your own turboshaft (s = 15: 1300 kW)

(a) T_rated = 594.0 N·m, limit 742.5 N·m, ultimate 1113.7 N·m. Ti-6Al-4V, 600 mm, 5° at limit:

| d_i/d_o | yield (mm) | ultimate (mm) | twist (mm) | chosen | kg/m | τ at limit |
|---|---|---|---|---|---|---|
| 0 | 19.52 | 21.51 | 33.28 | **33.28 (twist)** | 3.853 | 103 MPa |
| 0.8 | 23.27 | 25.64 | 37.96 | **37.96 (twist)** | 1.805 | 117 MPa |

(b) The ultimate criterion takes over when the twist requirement is relaxed beyond **28.6°** (solid) or **24.0°** (k = 0.8).

(c) Unless the drive-train dynamics genuinely allow more than about 25° of twist at limit torque, the stiffness requirement sets the shaft; it must therefore be specified explicitly, with its basis (a torsional natural-frequency requirement), and not left implicit behind a stress allowable.

**Check.** The diameters scale with power as P^(1/3) (strength) and P^(1/4) (twist): (1300/1400)^(1/4) × 33.90 = 33.28 mm.

**Common mistakes.** Sizing to the rated torque; forgetting the factor (1 − k⁴) in the twist formula.

## 5. Prandtl FD as the check

(a) 10 × 1 strip: FD **3.08822, 3.11447, 3.12105** at h = 0.1, 0.05, 0.025; series **3.12325**; bt³/3 = 3.33333 (the thin-strip formula without the end correction 1 − 0.630 b/a is 6.7 % high). Errors 3.50 × 10⁻², 8.78 × 10⁻³, 2.20 × 10⁻³: **order 2.0**.

(b) Hollow 2 × 1, wall 0.1: **0.214039, 0.215648, 0.216190, 0.216383** at h = t/2 … t/16. Successive differences 1.61 × 10⁻³, 5.42 × 10⁻⁴, 1.93 × 10⁻⁴ fall by 2.97 and 2.81: **order about 1.5**. The inner corners are re-entrant (270° of material around them); there the stress function's gradient (the shear stress) is singular, the solution is not smooth, and the 5-point scheme loses its second order near the corner.

(c) Bredt–Batho on the 1.9 × 0.9 mid-line rectangle: **0.208864**; ratio Bredt/FD = **0.965**. Consistent with the square tubes: a wall of 0.1 on a section of short side 1 sits between the lesson's t = 0.05 (0.9765) and t = 0.1 (0.9476) square cases.

**Check.** Extrapolating (b) with order 1.5 adds about 1.9 × 10⁻⁴ × 0.55 ≈ 1 × 10⁻⁴, so the converged J ≈ 0.2165, and the ratio moves by less than 0.001.

**Common mistakes.** Choosing h so that the wall is not a whole number of cells (the function refuses it); reading a converged result as exact when the order has dropped.

## 6. Torsion in practice

(a) An engine accelerating, a governor overshooting or a rotor suddenly loaded can briefly deliver more than its rated torque; the limit load must be the largest torque expected in service, and the ultimate factor 1.5 then covers the uncertainties beyond it. A real value comes from the engine's transient specification, from measured torque records, or from a drive-train dynamics model, and it can be well above 1.25 for a sudden stoppage.

(b) Torsional buckling of the thin wall (the wall of the k = 0.9 shaft is 2.2 mm on a 44 mm diameter), and whirling: the first bending critical speed of the shaft between its bearings must be well away from 20 900 rpm. A fatigue check at the splines is a third.

(c) In that bay the section is no longer closed: the shear flow cannot pass the hole, and the bay's torsional stiffness drops toward the open-section value, thousands of times lower (the lesson's 4886 for the whole box). The torque has to find another path, overloading the neighbouring bays and the ribs. A reinforced cut-out (a stiff frame or doubler around the hole, or a structural, fastened access panel that carries shear) restores the closed load path.

(d) Splines and keyways concentrate stress (a local stress several times the nominal one at a root radius), and a rotating shaft also sees alternating bending and torque ripple; fatigue life depends on the local alternating stress, so a crack starts at the root of a spline even when the body of the shaft is at 119 MPa.

## Check-yourself answers

1. ω = 628.3 rad/s; T = 100 000/628.3 = 159.2 N·m.
2. By symmetry: every diameter of a circle is equivalent, so nothing can make one part of the section warp differently from another. A square has corners and flats that are not equivalent, so its sections warp.
3. τ_max = 16T/(πd³); doubling d divides it by 8.
4. Area 1 − 0.8² = 36 %; J 1 − 0.8⁴ = 59 %.
5. T = 2A_m q; the shear flow is constant, so the stress q/t is highest where the wall is thinnest.
6. In the slit tube the shear flow cannot circulate round the cell and the torque is carried across the thickness with a lever arm of order t; the ratio is 3(r/t)² = 1200 for r/t = 20.
7. φ = 0 on the outer boundary; φ = constant on the boundary of each hole (with the constant fixed by the circulation condition).
8. Titanium (in the lesson, 1.873 against 2.425 kg/m at k = 0.8). For a twist-governed design the mass scales as ρ/√G for a fixed diameter ratio, so the ratio of density to the square root of the shear modulus decides it.
