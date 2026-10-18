# P10.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown.

## 1. The classic coarse lattice

**Assumptions.** Flat, untapered, 45° quarter-chord sweep, AR 5, one chordwise row.

(a) **CL_α = 3.4435 per radian** (0.0601 per degree), the value printed in the classic four-horseshoe-per-side textbook example.

(b)

| ny | uniform | Hough |
|---|---|---|
| 8 | 3.4435 | 3.1906 |
| 16 | 3.3045 | 3.1783 |
| 32 | 3.2354 | 3.1726 |
| 64 | 3.2020 | 3.1707 |

Observed orders from the last three: uniform ln(0.0690/0.0334)/ln 2 = **1.05**; Hough ln(0.00572/0.00194)/ln 2 = **1.56** (the swept wing's root kink, where the left and right bound legs meet at an angle, slows it below 2). Both head for about 3.17.

(c) The coarse answer is **8.6 %** above the converged 3.17. A coarse lattice is a teaching device; for design, refine until the change is below the accuracy you need and report the observed order.

**Common mistakes.** Comparing a coarse uniform lattice with a converged Hough one and calling the difference "the inset effect"; it is mostly resolution.

## 2. Derivations

(a) Vortex Γ at x = a, control point at x = a + d: tangency Γ/(2πd) = V∞α ⇒ Γ = 2πdV∞α; cl = 2Γ/(V∞c) = 4πdα/c. Equal to 2πα ⇒ **d = c/2**; the lift acts at the vortex, and thin-airfoil theory puts the flat plate's centre of pressure at c/4 ⇒ **a = c/4**.

(b) Panels of length c/2: vortices Γ₁ at x = c/8 and Γ₂ at 5c/8; control points at 3c/8 and 7c/8. A vortex upstream of a point induces downwash there, one downstream induces upwash. At 3c/8 (Γ₁ a distance c/4 upstream, Γ₂ c/4 downstream): (Γ₁ − Γ₂)/(2π c/4) = V∞α, so Γ₁ − Γ₂ = πcV∞α/2. At 7c/8 (Γ₁ 3c/4 upstream, Γ₂ c/4 upstream): Γ₁/(2π 3c/4) + Γ₂/(2π c/4) = V∞α, so Γ₁/3 + Γ₂ = πcV∞α/2. Adding the two equations: (4/3)Γ₁ = πcV∞α, so **Γ₁ = (3/4)πcV∞α, Γ₂ = (1/4)πcV∞α, Γ₁ + Γ₂ = πcV∞α and cl = 2πα** again, with the centre of pressure at (3/4)(c/8) + (1/4)(5c/8) = c/4.

(c) My(x₂) − My(x₁) = Σ(x₂ − x₁)F_z = (x₂ − x₁)ΣF_z. For each bound leg F = ρΓV × l with V = (cos α, 0, sin α), F_z = ρΓ cos α l_y, while the lift is ρΓ l_y; so ΣF_z = L cos α and ΔCm = CL cos α Δx/c_ref. Cm is therefore an affine function of the reference point and of CL, so dCm/dCL about x₂ = dCm/dCL about x₁ + cos α (x₂ − x₁)/c_ref (cos α ≈ 1 at small α): it vanishes at the single point **x_ac = x₁ − (dCm/dCL)₁ c_ref**.

(d) A 2-D vortex of strength s at (y_k, z_k) with its vector along +x induces v = (s/2πr²)(−(z − z_k), y − y_k) in the y–z plane. The far-field drag is the kinetic energy per unit length of the wake, D_i = ½ρ∫Γ w_T ds with w_T the *downwash normal to the strip* (positive down), i.e. −v·n; summed over strips, **D_i = −½ρ Σ G_j (v·n_j) Δs_j**.

**Common mistakes.** In (b), the sign of the upstream vortex's effect; in (c), forgetting cos α.

## 3. Forward sweep

(a)

| sweep | CL_α (/rad) | e | x_ac (m) | x_ac (MAC behind MAC LE) |
|---|---|---|---|---|
| −30° | 3.837 | 0.993 | −0.819 | 0.375 |
| 0° | 4.210 | 0.984 | 0.399 | 0.240 |
| +30° | 3.846 | 0.937 | 1.737 | 0.176 |

Forward and aft sweep cut the lift slope almost equally (cos Λ), but forward sweep keeps e high and moves the aerodynamic centre forward of the root.

(b) cl/CL is largest at **η = 0.02** (root) for −30° and 0° (peaks 1.334 and 1.186), and at **η = 0.42** for +30° (1.121). Aft sweep moves the stall outboard (towards the ailerons, the bad place); forward sweep keeps it at the root, which is one of its attractions.

(c) **Aeroelastic divergence**: when a forward-swept wing bends up under lift, the streamwise sections near the tip twist nose-up, which adds lift, which adds bending; above a critical dynamic pressure the loop diverges. Aft sweep has the opposite (stabilising) coupling. The lattice treats the wing as rigid; predicting divergence needs the structure's bending and torsion stiffness coupled to the aerodynamic loads (static aeroelasticity, a later phase).

**Common mistakes.** Measuring x_ac from the root LE and calling it a fraction of the MAC.

## 4. Your own convergence study (s = 15: n0 = 12)

(a) CL = **0.321268, 0.322644, 0.322979** for ny = 12, 24, 48; observed order ln(0.001376/0.000335)/ln 2 = **2.04**.

(b) Extrapolated CL = 0.322979 + 0.000335/3 = **0.323091**; CL_α = 0.323091/sin 4° = **4.632 per radian**, matching the lesson's 4.631 (nx 6, ny 40). The finest lattice is 1.1 × 10⁻⁴ (0.03 %) below it.

(c) With uniform spacing the values converge at order ≈ 0.9 from above and are several times further from the limit at the same ny (the lesson's table: 0.4099 at ny = 40 against 0.4036). With ny = n0 only the Hough lattice is trustworthy, because its error is small and its order is known.

**Common mistakes.** Extrapolating with an assumed order instead of the observed one; using α in degrees in sin.

## 5. The lattice against NACA TN-1669

(a)–(b) Per degree:

| AR | sweep | measured | VLM | VLM − measured | Helmbold–DATCOM |
|---|---|---|---|---|---|
| 5.16 | 0° | 0.0689 | 0.0698 | +0.0009 | 0.0751 |
| 5.16 | 45° | 0.0560 | 0.0561 | +0.0001 | 0.0592 |
| 5.16 | 60° | 0.0434 | 0.0430 | −0.0004 | 0.0452 |
| 2.61 | 0° | 0.0486 | 0.0508 | +0.0022 | 0.0541 |
| 2.61 | 45° | 0.0463 | 0.0448 | −0.0015 | 0.0462 |
| 2.61 | 60° | 0.0360 | 0.0367 | +0.0007 | 0.0377 |
| 1.34 | 0° | 0.0280 | 0.0323 | +0.0043 | 0.0333 |
| 1.34 | 45° | 0.0315 | 0.0310 | −0.0005 | 0.0309 |
| 1.34 | 60° | 0.0286 | 0.0282 | −0.0004 | 0.0275 |

Within the ±0.0005 reading uncertainty: AR 5.16 at 45° and 60°, AR 1.34 at 45° and 60°; just outside at AR 5.16 unswept (+1.3 %, plausibly the section slope below 2π) and AR 2.61 at 60°. Clearly outside: AR 2.61 unswept and 45° (±3–5 %) and AR 1.34 unswept (+15 %).

(c) A low-aspect-ratio rectangular wing sheds strong vortices from its *side edges*, which separate and roll up over the upper surface; their suction adds lift that grows nonlinearly with α (∝ sin²α). The measured lift curve is then curved, and a slope read at small CL can differ from the linear theory's; moreover the section slope at R ≈ 1.6 × 10⁶ is below 2π. The VLM, linear and attached, has neither effect. (Polhamus' leading-edge suction analogy is the classical way to add it.)

**Common mistakes.** Comparing per-radian VLM slopes with per-degree data; treating digitised data as exact.

## 6. Lattice or line?

(a) The lattice: it resolves the chordwise extent of the tips and converges with a known order. At cruise CDi = 0.00431 (95.7 N) against 0.00448 (99.6 N) from the lifting line: **about 4 N** of 783 N total, smaller than the uncertainty in the aircraft's Oswald factor (`P10.09`), so either is acceptable for performance; the lattice is the better number.

(b) A planar wing cannot have e > 1. Uniform spacing puts the tip vortex exactly at the tip; on a coarse lattice this misplaces the tip's load and the Trefftz sum overshoots (the lesson's table: e = 1.070 at ny = 10). Use Hough's quarter-panel inset (or a much finer lattice) and report the convergence.

(c) (1) Profile drag and the drag polar: `P10.06`, `P10.13`. (2) cl_max, stall progression and post-stall behaviour: `P10.15` (with `P10.12` for where it starts). (3) Surface pressures with thickness, and the section's own pitching moment: `P10.05`, `P10.06`. Also compressibility, aeroelastic deformation and propeller slipstream.

## Check-yourself answers

1. The bound vortex at each panel's quarter chord, the control point at its three-quarter chord; it reproduces the 2-D flat plate's lift slope 2π and centre of pressure at c/4.
2. The normal velocity at control point i induced by a horseshoe j of unit strength: AIC_ij = v_j(P_i)·n_i.
3. The near-field (bound-leg) induced drag needs the induced velocity at the bound legs, which converges slowly and is sensitive to the lattice; the Trefftz-plane drag depends only on the shed circulation and converges smoothly.
4. The point about which the pitching moment does not change with lift: x_ac = x_ref − (dCm/dCL) c_ref.
5. It concentrates the tip's fall-off of load into the last strip, whose error is first order in the strip width; Hough's quarter-panel inset fixes the placement and gives order 2.
6. From three solutions on lattices refined by a constant ratio r: p = ln[(f₁ − f₂)/(f₂ − f₃)]/ln r.
7. The lifting line ignores that the load falls off over a chord-length near each tip; that region is a fraction 1/AR of the span, so the difference is about −0.43/AR.
8. Any two of: stall and separation; viscous (profile) drag; thickness effects on the pressures; separated leading-edge or side-edge vortex lift; compressibility; wake roll-up.
