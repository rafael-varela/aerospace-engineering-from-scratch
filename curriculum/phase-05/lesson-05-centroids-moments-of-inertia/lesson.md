# P05.05  Centroids & Moments of Inertia

> A beam does not know its material's strength until it knows its shape: area, centroid and second moments are the numbers that turn a cross-section into a stiffness.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the centroid of a composite area and of a set of masses (an aircraft's centre of gravity in % MAC) as a weighted average. 2. Compute the second moments Ix, Iy and the product Ixy of a composite section with the parallel-axis theorem, including holes as negative parts. 3. Derive and implement the exact area integrals of any polygon from Green's theorem, and use them as an independent check of the composite method. 4. Rotate second moments to other axes, find the principal second moments and axes of an unsymmetric section by Mohr's circle, and verify them as the eigenvalues of a 2 × 2 matrix. 5. Quantify the error of thin-walled and boom idealisations of stringers and wing boxes, and judge when they are good enough. |
| **Prerequisites** | `P05.01` (resultant and centroid of a distributed load), `P01.10` (multiple integrals), `P01.12` (Green's theorem), `P01.16` (eigenvalues of a symmetric matrix), `P01.15` (rotation matrices). Skill check: you can evaluate ∫∫ y² dA over a rectangle by hand and compute the eigenvalues of a 2 × 2 symmetric matrix with `numpy.linalg.eigh`. |
| **Engineering question** | A light-aircraft wing spar is an extruded I-section, its skin is stiffened by Z-section stringers, and the wing box is often analysed as a set of point areas. Where are the centroids, what are the second moments and principal axes of these sections, how much does the I-section gain over a solid bar of the same height, and how wrong are the thin-walled and boom idealisations that structural analysts use every day? |
| **Assumptions and limits** | Sections made of straight-sided pieces (curves are approximated by polygons, with the error measured); sharp corners (no fillet radii); nominal thicknesses; one homogeneous material per section unless stated (the weight-and-balance part uses masses directly). Section properties are geometry only: how they enter stress and deflection is `P05.10` to `P05.12`. |
| **Analytical work** | Centroid as a weighted average; first and second moments of area; the parallel-axis theorem for I and for Ixy; composite tables with negative parts; Green's theorem for the polygon integrals; rotation of axes and its invariants; principal axes by Mohr's circle and by eigenvalues; radius of gyration and section modulus; thin-walled approximations and their (t/b)² error; the boom idealisation of a stiffened box; centre of gravity in % MAC. |
| **Computational deliverable** | `starter/main.py`: `polygon_integrals`, `section_from_integrals`, `polygon_section`, `rect`, `circle`, `boom`, `composite_section`, `rotate_axes`, `principal_axes`, `principal_axes_eig`, `radius_of_gyration`, `section_moduli`, `z_section_thin`, `wing_box_parts`, `composite_centroid`, `cg_percent_mac` (the `Section` and `Part` dataclasses, the shape builders `rect_polygon`, `circle_polygon`, `i_section_*`, `z_section_polygon`, `z_section_parts`, `angle_polygon`, and `main` are given). Tests prove: the rectangle and triangle closed forms; vertex order, start and translation do not change the polygon results; the parallel-axis theorem on the raw integrals; inscribed polygons converge to the circle at second order; composite and polygon agree for the I-spar and the Z-stringer (and the I-spar matches its textbook formula); holes by negative parts and negative polygons agree; rotation preserves the invariants and a quarter turn swaps axes; rotated second moments equal those of the rotated polygon; principal axes agree with `eigh` and zero the product; symmetric sections have principal axes along x and y; the thin-walled error falls like t²; the wing-box idealisations behave as expected; section moduli, radius of gyration, centroid of masses and % MAC. |
| **Practical activity** | **computational experiment**: the error of the thin-walled Z-section formulas against wall thickness, the principal angle of a Z against flange width, the convergence of a polygon to a circle, and three idealisations of a wing box. |
| **Equipment** | None. A computer with the environment from `P00.01`. A sheet of card, scissors and a pin are a pleasant way to find a centroid by suspension, but no data from it are used. |
| **Safety** | None beyond normal computer use (and care with scissors if you try the card demonstration). |
| **Evidence** | Completed `main.py`, passing tests, `section_properties.png`, a hand composite table for the I-spar, the derivation of the polygon area and centroid formulas, the four experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: a correct composite table by hand; both routes (composite and polygon) agree in code; principal axes found and checked by eigenvalues; the thin-walled and boom errors stated as numbers with their trend. |
| **References** | `@libretexts-baker-haynes-engineering-statics` Ch. 7 "Centroids and Centers of Gravity": 7.1 Weighted Averages, 7.2 Center of Gravity, 7.4 Centroids, 7.5 Centroids using Composite Parts, 7.7 Centroids using Integration (required-free); Ch. 10 "Moments of Inertia": 10.1 Integral Properties of Shapes, 10.2 Moments of Inertia of Common Shapes, 10.3 Parallel Axis Theorem, 10.4 Moment of Inertia of Composite Shapes, 10.6 Radius of Gyration, 10.7 Products of Inertia (required-free). `@mit-ocw-16-001-unified-materials-structures` lecture notes "Moment of Intertia [sic] of Beam Cross Section, part 01" and "part 02" (optional-free: the aerospace view of the same quantities). |
| **Extensions** | ★ Add fillet radii to the I-spar polygon (each corner a quarter-circle polygon) and measure their effect on A and Ix. ★ Compute the polar moment J = Ix + Iy and show that it is the same about every pair of axes through the centroid. ★ Extend `composite_section` to parts of different materials by weighting each with E_i/E_ref (a transformed section, used in `P05.10`). |
| **Time** | unmeasured |

## 1. The engineering problem

When a wing bends, the material far from the neutral axis does most of the work. A wing spar is therefore made as an **I-section**: two flanges (spar caps) far apart, joined by a thin web. A wing skin is stiffened by **stringers**, often Z-sections, which are unsymmetric and bend sideways when loaded vertically. And a whole wing box is often analysed as a handful of **booms**, point areas at the spar caps and stringers, with the skin either ignored or lumped into the booms.

All of this rests on a few area integrals of the cross-section: the area A, the **centroid** (through which the neutral axis passes in pure bending, `P05.10`), the **second moments** Ix = ∫y² dA and Iy = ∫x² dA, and the **product** Ixy = ∫xy dA. Their mass counterpart, the centroid of a set of masses, is the aircraft's **centre of gravity**, which must lie inside limits stated in % of the mean aerodynamic chord.

This lesson builds a **section property calculator** two independent ways, by composite parts with the parallel-axis theorem and by exact polygon integrals from Green's theorem, and uses it on a 150 mm × 50 mm extruded I-spar, a 25 mm Z-stringer, an unequal angle, a thin tube and a 400 mm × 120 mm stiffened wing box. Headline results: the I-spar's flanges carry **82.6 %** of its Ix with **59.2 %** of its area, and a solid bar of the same width would need **4.8 times the area** for the same Ix; the Z-stringer's principal axes are turned **28.5°** from its web; and the boom idealisation without skin under-estimates the wing box's Ix by **45 %**, while lumping the skin into the booms brings it within **1.3 %**.

Success criteria, written first: composite and polygon agree to round-off; principal axes from Mohr's circle equal the eigenvectors of the second-moment matrix; every idealisation has its error measured against an exact calculation.

## 2. Understand and derive

### 2.1 Centroid: a weighted average

The centroid of an area is the average position of its area (Baker & Haynes 7.1, 7.4):

x̄ = (1/A) ∫ x dA = Q_y/A,   ȳ = (1/A) ∫ y dA = Q_x/A,   (1)

where Q_x = ∫ y dA and Q_y = ∫ x dA are the **first moments** of the area. For a body made of parts with known areas and centroids, the integrals split into sums (7.5): x̄ = Σ A_i x̄_i / Σ A_i, with a hole entered as a negative area. The same formula with masses instead of areas gives the **centre of gravity** of a set of masses (7.2), the basis of every aircraft's weight-and-balance sheet.

An illustrative light aircraft: empty mass 700 kg at 2.10 m aft of the datum, pilot 80 kg and passenger 75 kg at 2.00 m, 100 L of fuel (72 kg) at 2.25 m and 20 kg of baggage at 3.00 m. Then m = 947 kg and x_cg = 2002.0/947 = **2.1140 m**. With the leading edge of the mean aerodynamic chord (LEMAC) at 1.80 m and a MAC of 1.50 m, the CG is at (2.1140 − 1.80)/1.50 = **20.94 % MAC**. (These numbers are made up for the example; a real aircraft's limits come from its flight manual.)

### 2.2 Second moments and the parallel-axis theorem

The **second moment of area** (often called the area moment of inertia) about the x-axis, Ix = ∫ y² dA, measures how far the area is spread from that axis; bending stiffness about x is E·Ix (`P05.10`). For a rectangle b wide and h high about its centroidal axis, Ix = bh³/12 (10.2): the cube is why depth is cheap stiffness.

About an axis parallel to a centroidal one and a distance d from it, write y = y_c + d and expand: ∫(y_c + d)² dA = ∫y_c² dA + 2d ∫y_c dA + d²A. The middle integral is zero because y_c is measured from the centroid. So (10.3)

**I = I_c + A d²**,  and likewise  **Ixy = Ixy_c + A d_x d_y**.   (2)

The theorem only works *from the centroid*: moving from one non-centroidal axis to another needs two steps. The product Ixy is zero for any section symmetric about x or y (each element x y dA has a mirror image −x y dA).

**The I-spar by hand** (150 mm high, flanges 50 × 6 mm, web 3 mm thick; centroid at the centre by symmetry):

| part | A (mm²) | d (mm) | own I_c (mm⁴) | A d² (mm⁴) | Ix (mm⁴) |
|---|---|---|---|---|---|
| top flange 50 × 6 | 300 | 72 | 900 | 1 555 200 | 1 556 100 |
| bottom flange 50 × 6 | 300 | −72 | 900 | 1 555 200 | 1 556 100 |
| web 3 × 138 | 414 | 0 | 657 018 | 0 | 657 018 |
| **total** | **1 014** | | | | **3 769 218** |

The flanges' own second moments (900 mm⁴ each) are negligible: almost all their contribution is A d². This is why aircraft structures are idealised as booms. The textbook formula for a symmetric I, (bh³ − (b − t_w)(h − 2t_f)³)/12, gives the same 3 769 218 mm⁴. The flanges carry **82.6 %** of Ix with **59.2 %** of the area. A solid bar 50 mm wide with the same Ix would be 96.7 mm high and have **4 836 mm²**, 4.8 times the area; the full 150 × 50 mm bar would have Ix = 14 062 500 mm⁴ but weigh 7.4 times as much.

Two derived quantities complete the picture (10.6): the **radius of gyration** r = √(I/A), the distance at which all the area would give the same I (r_x = 60.97 mm, r_y = 11.12 mm for the spar: it will buckle about y, `P05.15`), and the **section modulus** S = I/c, with c the distance to the extreme fibre (S = 50 256 mm³), which turns a bending moment into the largest stress in `P05.10`. In 7075-T6 aluminium, with a typical density of 2 810 kg/m³ (handbooks such as MMPDS, formerly MIL-HDBK-5, give about 0.10 lb/in³, 2 800 kg/m³), the spar weighs **2.849 kg/m**.

### 2.3 Exact integrals of a polygon: Green's theorem

Any section with straight edges is a polygon, and its area integrals have exact closed forms. Green's theorem (`P01.12`) turns an area integral into a boundary integral: for example A = ∮ x dy = ½ ∮ (x dy − y dx). Along a straight edge from (x_i, y_i) to (x_{i+1}, y_{i+1}), with c_i = x_i y_{i+1} − x_{i+1} y_i, each integral becomes one term per edge:

A = ½ Σ c_i,   Q_y = ⅙ Σ (x_i + x_{i+1}) c_i,   Q_x = ⅙ Σ (y_i + y_{i+1}) c_i,

I_yy = (1/12) Σ (x_i² + x_i x_{i+1} + x_{i+1}²) c_i,   I_xx = (1/12) Σ (y_i² + y_i y_{i+1} + y_{i+1}²) c_i,

I_xy = (1/24) Σ (x_i y_{i+1} + 2x_i y_i + 2x_{i+1} y_{i+1} + x_{i+1} y_i) c_i.   (3)

These are about the **origin**; the centroidal values follow from (1) and (2) used backwards: Ix = I_xx − A ȳ², Iy = I_yy − A x̄², Ixy = I_xy − A x̄ ȳ. The sums are signed: an anticlockwise polygon gives A > 0, a clockwise one A < 0, so the code normalises the sign. A hole is a polygon subtracted. The formulas are exact for straight edges: the I-spar as one 12-vertex polygon gives exactly the composite table's numbers, an independent check with no shared code. A curve is approximated by an inscribed polygon, with an error that falls as 1/n² (Section 4).

### 2.4 Rotated axes and principal axes

For axes u, v turned anticlockwise by θ from x, y through the centroid, substituting u = x cos θ + y sin θ, v = −x sin θ + y cos θ into the integrals gives

I_u = Ix cos²θ − 2 Ixy sin θ cos θ + Iy sin²θ,
I_v = Ix sin²θ + 2 Ixy sin θ cos θ + Iy cos²θ,
I_uv = (Ix − Iy) sin θ cos θ + Ixy (cos²θ − sin²θ).   (4)

Two **invariants** follow: I_u + I_v = Ix + Iy (the polar moment) and I_u I_v − I_uv² = Ix Iy − Ixy². Equations (4) are a rotation of a symmetric tensor, the same algebra as stress transformation (`P05.13`), and they plot as **Mohr's circle**: centre (Ix + Iy)/2, radius R = √(((Ix − Iy)/2)² + Ixy²). The **principal second moments** are the extremes, **I_max,min = (Ix + Iy)/2 ± R**, on the **principal axes** where I_uv = 0:

**tan 2θ_p = −2 Ixy/(Ix − Iy).**   (5)

A section with an axis of symmetry has its principal axes along and across it. An unsymmetric section such as a Z or an angle does not, and that matters: bent about a non-principal axis, it deflects sideways as well (unsymmetric bending, `P05.10`); and it buckles about its axis of least I (`P05.15`).

**An independent check by eigenvalues.** The second moment about an axis with unit normal n = (n_x, n_y) is I = ∫(n·r)² dA = nᵀ C n with **C = [[Iy, Ixy], [Ixy, Ix]]** (careful: x² dA goes with n_x²). The principal values are the eigenvalues of C, and its eigenvectors are the principal axes' *normals* (`P01.16`). `numpy.linalg.eigh` gives them with no trigonometry.

**The Z-stringer** (web 25 mm, flanges 15 mm, all 1.6 mm thick): A = 82.88 mm², Ix = 7 962.3, Iy = 3 062.2, **Ixy = 3 762.7 mm⁴** (both flanges lie in quadrants where xy > 0). Mohr's circle: centre 5 512.2, R = 4 490.1, **I_max = 10 002.3 and I_min = 1 022.1 mm⁴** at **θ_p = −28.465°**; `eigh` gives the same values, and rotating by θ_p makes I_uv = 4 × 10⁻¹³ mm⁴. The minimum second moment is only an eighth of Ix. **The angle** 30 × 20 × 2 mm (heel at the origin): centroid (4.75, 9.75) mm, Ix = 8 922, Iy = 3 242, Ixy = −3 150 mm⁴ (the legs lie where xy < 0 relative to the centroid), principal 10 323.2 and 1 840.8 mm⁴ at +23.981°.

### 2.5 Thin-walled sections and booms

Aircraft sections are thin: t/b of 0.05 to 0.15. Treating the wall as a line along its mid-line and dropping terms in t³ gives short formulas. For the Z, with web length d = h − t and flange length f = b − t/2 measured on the mid-line:

Ix ≈ t d³/12 + 2 f t (d/2)²,   Iy ≈ 2 t f³/3,   Ixy ≈ 2 f t (f/2)(d/2).   (6)

For the 1.6 mm Z the errors against the exact polygon are −0.42 %, −0.26 % and +0.32 %. The dropped terms are of relative order (t/b)², so the error falls like t²: halving the wall quarters it.

The **boom idealisation** goes further: the spar caps and stringers become point areas (their own I_c dropped, which §2.2 showed is negligible), and the thin skin and webs are either dropped, or their bending effectiveness is lumped into the booms. For a skin panel of width w and thickness t between two booms whose direct stresses are σ₁ and σ₂, the classic lumping adds t w/6 × (2 + σ₂/σ₁) to boom 1 (σ₂/σ₁ = 1 for two booms at the same height in bending about x; −1 for the two ends of a vertical web). The rule makes a boom pair have the same second moment as the panel between them. For the 400 × 120 mm box (skins 1.6 mm, webs 2 mm, four 300 mm² caps, four 60 mm² stringers per skin): full Ix = 11 087 781 mm⁴; **booms alone 6 048 000 mm⁴ (−45.5 %)**; **skin lumped into booms 11 232 000 mm⁴ (+1.30 %)**. The skin is not negligible in bending, but it is easy to account for. The box weighs 9.55 kg/m in 2024-T3 (typical density 2 780 kg/m³; MMPDS gives about 0.10 lb/in³).

### 2.6 Where the model stops being valid

- **Fillet radii and tolerances.** Extrusions have fillets at the web–flange junction, and sheet thickness has a tolerance. A ±0.1 mm change on the 1.6 mm Z changes Ix by about ±5 %, ten times the thin-walled error.
- **Holes and cut-outs.** Rivet and fastener holes reduce the section at the hole station; the net section governs the stress there but not the overall stiffness.
- **Effective width.** Thin skins in compression buckle between stringers and then carry less than their full width (`P05.15` and later phases); the lumped-boom result assumes the skin is fully effective.
- **Different materials.** Composite parts of different stiffness need transformed sections (Extensions).

## 3. Build the computational model

Open `starter/main.py`. The `Section` dataclass (A, cx, cy, Ix, Iy, Ixy, all centroidal, SI units) and the `Part` dataclass (a piece's own centroidal properties and a sign, −1 for a hole) are given, with the shape builders and `main`.

### Step 1: polygon integrals

`polygon_integrals(xy)` vectorises (3) with `np.roll` and returns (A, Q_x, Q_y, I_xx, I_yy, I_xy) about the origin, with the sign normalised so that A > 0. `section_from_integrals` applies (1) and (2) backwards; `polygon_section(outer, holes)` subtracts the holes' integrals before that step.

### Step 2: composite sections

`rect`, `circle` and `boom` return `Part`s; `composite_section(parts)` is (1) and (2):

```python
A = np.array([p.A for p in parts]) * s
cx = np.sum(A * [p.cx for p in parts]) / A.sum()
Ix = sum(p.sign * (p.Ix + p.A * (p.cy - cy) ** 2) for p in parts)
Ixy = sum(p.sign * (p.Ixy + p.A * (p.cx - cx) * (p.cy - cy)) for p in parts)
```

### Step 3: rotation and principal axes

`rotate_axes` is (4); `principal_axes` returns (I_max, I_min, θ_p) from Mohr's circle, with θ_p = ½ atan2(−2Ixy, Ix − Iy) so that the axis returned is the one of I_max; `principal_axes_eig` returns `eigh(C)`. `radius_of_gyration` and `section_moduli` complete the set.

### Step 4: idealisations and masses

`z_section_thin` is (6); `wing_box_parts` builds the box with `idealisation` = "full", "booms" or "lumped"; `composite_centroid` is (1) for masses, and `cg_percent_mac` converts to % MAC.

Running `main()` prints every number of this lesson and writes `section_properties.png` (the sections with their centroids and principal axes; the thin-walled error against t; the polygon convergence to a circle).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-05-centroids-moments-of-inertia`. The tests check:

1. **Closed forms**: rectangle bh³/12, hb³/12 and Ixy = 0; right triangle bh³/36 and Ixy = −b²h²/72.
2. **Invariance**: reversing the vertex order, starting at another vertex or translating the polygon leaves the centroidal properties unchanged.
3. **Parallel-axis theorem** on the raw integrals of an offset rectangle.
4. **Convergence**: inscribed n-gons converge to π r⁴/4 at second order.
5. **Two routes agree**: composite = polygon for the I-spar and the Z (and the I-spar matches the textbook formula); a hole as a negative part = a hole as a subtracted polygon.
6. **Rotation**: the invariants hold at every angle; a quarter turn swaps Ix and Iy and changes the sign of Ixy; (4) equals the integrals of the rotated polygon.
7. **Principal axes**: Mohr = `eigh`; I_uv = 0 on the principal axes; symmetric sections give θ_p = 0 or 90°.
8. **Thin-walled error** falls like t² (the error ratio is 4 when t doubles, within 5 %, for very thin walls).
9. **Wing box**: booms only < full; lumped within 5 % of full and booms only more than 5 % off; an unknown idealisation is refused.
10. **Moduli, gyration, centroid of masses, % MAC.**
11. `main` writes the figure.

### Compare

The composite and polygon routes share no code beyond the `Section` container, and they agree to 10⁻¹² relative for every straight-sided section in the lesson. Mohr's circle and `eigh` agree to round-off. The thin-walled and boom idealisations do not agree with the exact values, and are not supposed to: their errors are the results of the experiment.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are computations; the sections are nominal geometry. Write each prediction down before running it.

1. **Thin-walled error.** For the Z with h = 25 mm and b = 15 mm, compute the error of (6) for t = 0.5, 1.0, 1.6, 2.5 and 4.0 mm. Predict how it scales with t.
2. **Principal angle.** For the Z with h = 25 mm and t = 1.6 mm, vary the flange width b from 5 to 30 mm and record θ_p and I_min/I_max. Predict the trend of θ_p.
3. **Polygon convergence.** Approximate a 20 mm × 1 mm tube by inscribed n-gons (outer and inner), n = 8 to 256, and compare Ix with the exact π(D⁴ − d⁴)/64. Predict the order.
4. **Wing-box idealisations.** Compare full, booms-only and lumped Ix for the box of §2.5.

Raw data format: CSV tables `t_mm, t_over_b, err_Ix_pct, err_Iy_pct, err_Ixy_pct` (1); `b_mm, theta_p_deg, Imin_over_Imax` (2); `n, Ix_mm4, rel_error` (3); `model, Ix_mm4, error_pct` (4); plus `section_properties.png`.

What you should find (model results): (1) Ix error −0.038 % at 0.5 mm, −0.158 % at 1.0, −0.423 % at 1.6, −1.105 % at 2.5 and −3.159 % at 4.0 mm: doubling t multiplies the error by almost exactly four, the (t/b)² law. Iy and Ixy behave the same way with smaller coefficients. (2) θ_p = −5.17°, −16.27°, −28.47°, −39.66°, −48.61°, −55.33° for b = 5 to 30 mm: the wider the flanges, the further the principal axes turn from the web. I_min/I_max peaks near b = 20 mm (0.116) and is never above 0.12: a Z is always much weaker about its minor axis. (3) Relative error 1.876 × 10⁻¹ (8-gon), 5.02 × 10⁻², 1.28 × 10⁻², 3.21 × 10⁻³, 8.03 × 10⁻⁴, 2.008 × 10⁻⁴ (256-gon): observed order 1.90, 1.98, 1.99, 2.00, 2.00. An inscribed polygon always under-estimates (it lies inside the circle). (4) Booms only −45.5 %, lumped +1.30 %.

## 5. Validate and improve

| quantity | hand / closed form | code | difference | cause |
|---|---|---|---|---|
| I-spar Ix (mm⁴) | 3 769 218 (table; textbook formula) | 3 769 218.0 (composite and polygon) | 10⁻¹² | round-off |
| Z principal values (mm⁴) | Mohr: 10 002.3, 1 022.1 | `eigh`: 10 002.34, 1 022.15 | 10⁻¹² | round-off |
| tube Ix (mm⁴) | 2 700.984 | 2 700.442 (256-gon) | 2 × 10⁻⁴ | polygon approximation, O(1/n²) |
| Z Ix (mm⁴) | thin-walled 7 928.7 | exact 7 962.3 | −0.42 % | t³ terms dropped |
| box Ix (mm⁴) | lumped 11 232 000 | full 11 087 781 | +1.30 % | lumping rule is approximate |
| CG (% MAC) | 2002.0/947 → 20.94 | 20.94 | 0 | — |

The calculator agrees with every closed form; the idealisations' errors are measured and small enough for preliminary design except the booms-only model. The next improvement would be geometric realism: fillet radii, actual thickness from a measured part, and the net section at fastener holes.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand composite table for the I-spar; the derivation of A and x̄ for a polygon from Green's theorem; `section_properties.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the section is turned upside down?** Ix and Iy are unchanged, Ixy changes sign for a reflection (not for a 180° rotation). A sign error in Ixy flips the principal axis to the wrong side of the web: always check against `eigh` or a sketch of where the area lies.
- **What if the thickness goes to zero?** A, Ix, Iy and Ixy all go to zero in proportion to t, and the thin-walled error vanishes as t²; the principal angle tends to a finite limit set by the mid-line geometry.
- **What if a part is entered with the wrong sign (a hole as material)?** The area grows instead of shrinking; a check of the net area against a sketch catches it at once.
- **Which assumption fails first in the real world?** Nominal thickness: the tolerance on sheet and extrusion thickness gives a larger change in Ix than the thin-walled approximation.
- **What measurement could falsify the model?** Weigh a measured length of extrusion and compare mass/(ρL) with A; or bend it and back out EI (`P05.18`).
- **Does the result survive a change of method?** Yes: composite and polygon routes agree to round-off, and Mohr's circle and the eigenvalue route give the same principal axes.

## Check yourself

1. Define the centroid of an area as a weighted average and write it for a set of parts with one hole.
2. State the parallel-axis theorem for I and for Ixy, and explain why it only works from the centroid.
3. Why do the flanges of an I-section carry most of its Ix, and why can their own centroidal second moments be neglected?
4. Write the shoelace (Green's theorem) formula for the area of a polygon. What happens to its sign if the vertices are listed clockwise?
5. Why is Ixy zero for a section with an axis of symmetry, and what does a non-zero Ixy tell you?
6. Give the principal second moments of a section from Ix, Iy, Ixy, and say how to check them with an eigenvalue routine.
7. A thin-walled formula is 0.4 % in error at t = 1.6 mm. Roughly what will it be at t = 3.2 mm, and why?
8. An aircraft of 947 kg has its CG at 2.114 m; 20 kg of baggage is added at 3.00 m. Where is the new CG?
