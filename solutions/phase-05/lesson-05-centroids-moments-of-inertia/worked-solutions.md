# P05.05  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` (and the probe script that calls it) and given to the figures shown. Dimensions are nominal; material densities are typical values, not measurements.

## 1. A T-section by hand (s = 15: h = 55 mm)

**Assumptions.** Sharp corners, nominal dimensions; origin at the bottom of the web.

(a)

| part | A (mm²) | ȳ (mm) | Aȳ (mm³) | I_c (mm⁴) | d = ȳ − Ȳ (mm) | A d² (mm⁴) |
|---|---|---|---|---|---|---|
| flange 60 × 8 | 480 | 51.0 | 24 480 | 2 560 | 10.177 | 49 715.9 |
| web 6 × 47 | 282 | 23.5 | 6 627 | 51 911.5 | −17.323 | 84 622.7 |
| **total** | **762** | | **31 107** | | | |

Ȳ = 31 107/762 = **40.823 mm** above the base. Ix = 2 560 + 49 715.9 + 51 911.5 + 84 622.7 = **188 810.1 mm⁴**.

(b) Both parts are centred on the y-axis: Iy = 8 × 60³/12 + 47 × 6³/12 = 144 000 + 846 = **144 846.0 mm⁴**. The section is symmetric about y, so every element x y dA has a mirror −x y dA: **Ixy = 0**.

(c) c_top = 55 − 40.823 = 14.177 mm, c_bottom = 40.823 mm. **S_top = 13 317.9 mm³**, **S_bottom = 4 625.1 mm³**; **r_x = √(188 810.1/762) = 15.74 mm**. The bottom fibre (the tip of the web), farther from the centroid, sees the larger bending stress, 2.9 times the top's.

(d) `composite_section` and `polygon_section` on the outline (−3, 0), (3, 0), (3, 47), (30, 47), (30, 55), (−30, 55), (−30, 47), (−3, 47) mm both give A = 762.000 mm², Ȳ = 40.8228 mm, Ix = 188 810.1 mm⁴, Iy = 144 846.0 mm⁴.

**Check.** The centroid lies between the parts' centroids (23.5 and 51 mm), nearer the larger area. Dimensions: mm⁴ for I, mm³ for S.

**Common mistakes.** Adding the parts' I_c about their own centroids without A d²; measuring d from the base instead of from the section centroid; double-counting the overlap of flange and web (here the web is taken as 47 mm, up to the flange).

## 2. Derivations

(a) On the edge, x = x_i + τΔx, y = y_i + τΔy, τ ∈ [0, 1], so x dy − y dx = (x_i + τΔx)Δy dτ − (y_i + τΔy)Δx dτ = (x_iΔy − y_iΔx) dτ, and x_iΔy − y_iΔx = x_i y_{i+1} − x_{i+1} y_i = c_i. Integrating over τ and summing, **A = ½ Σ c_i** (the shoelace formula). **Q_y**: the triangle with vertices at the origin, vertex i and vertex i + 1 has signed area c_i/2 (the same cross-product) and its centroid's x is (0 + x_i + x_{i+1})/3. The polygon is the signed sum of these triangles (areas outside it cancel between edges), so Q_y = Σ (c_i/2)(x_i + x_{i+1})/3 = **⅙ Σ (x_i + x_{i+1}) c_i**. Check: a unit square (0,0), (1,0), (1,1), (0,1) gives c = 0, 1, 1, 0, A = 1 and Q_y = ⅙(2 + 1) = ½, so x̄ = ½.

(b) With y = y_c + d, where y_c is measured from the centroidal axis: Ix = ∫(y_c + d)² dA = ∫y_c² dA + 2d ∫y_c dA + d² A. The middle integral is the first moment about the centroidal axis, zero by definition of the centroid, so **Ix = Ix_c + A d²**. Likewise ∫(x_c + d_x)(y_c + d_y) dA = Ixy_c + d_x ∫y_c dA + d_y ∫x_c dA + A d_x d_y = **Ixy_c + A d_x d_y**. If the starting axis is not centroidal, the middle term does not vanish: the theorem gives a wrong answer unless one first moves back to the centroid.

(c) With u = x cos θ + y sin θ and v = −x sin θ + y cos θ (anticlockwise axes), I_u = ∫v² dA = Ix cos²θ − 2Ixy sin θ cos θ + Iy sin²θ, I_v = ∫u² dA = Ix sin²θ + 2Ixy sin θ cos θ + Iy cos²θ, I_uv = ∫uv dA = (Ix − Iy) sin θ cos θ + Ixy(cos²θ − sin²θ). Adding gives I_u + I_v = Ix + Iy. With double angles, I_u = c + a cos 2θ − Ixy sin 2θ, I_v = c − a cos 2θ + Ixy sin 2θ, I_uv = a sin 2θ + Ixy cos 2θ, where c = (Ix + Iy)/2 and a = (Ix − Iy)/2; then I_u I_v − I_uv² = c² − (a cos 2θ − Ixy sin 2θ)² − (a sin 2θ + Ixy cos 2θ)² = c² − a² − Ixy² = Ix Iy − Ixy². **I_uv = 0** when **tan 2θ_p = −2Ixy/(Ix − Iy)**. For an axis through the centroid with unit normal n, a point r is at distance n·r = n_x x + n_y y from it, so I = ∫(n_x x + n_y y)² dA = n_x² Iy + 2n_x n_y Ixy + n_y² Ix = **nᵀCn** with C = [[Iy, Ixy], [Ixy, Ix]]. Its stationary values on |n| = 1 are the eigenvalues of C.

**Common mistakes.** In (c), putting Ix where Iy belongs in C (x² dA multiplies n_x²); in (a), forgetting that the area sum is signed.

## 3. A helicopter blade D-spar (s = 15: rod diameter 9 mm)

**Assumptions.** Nominal geometry, the hole is the inner circle of radius 27.5 mm cut off by the web's inner face at x = −2 mm; sharp corners; one material for the spar.

(a) With 2000 points per arc: **A = 335.70 mm²**, centroid **x̄ = −12.646 mm** (ȳ = 0 by symmetry), **Ix = 121 152 mm⁴**, **Iy = 39 957.5 mm⁴**, Ixy = 0.

| points per arc | A (mm²) | x̄ (mm) | Ix (mm⁴) |
|---|---|---|---|
| 50 | 335.44 | −12.639 | 120 970 |
| 200 | 335.69 | −12.645 | 121 141 |
| 2000 | 335.70 | −12.646 | 121 152 |

The 200-point result is within 0.01 % of the 2000-point one.

(b) Mid-line radius r_m = 28.75 mm: arc area π r_m t = 225.8 mm² at x = −2r_m/π = −18.30 mm; web area 2r_m × 2 = 115.0 mm² at x = −1 mm. Total **A ≈ 340.8 mm²** (+1.5 %), **x̄ ≈ −12.46 mm** (0.18 mm aft of the exact). The thin-walled Ix, t r_m³ π/2 + 2(2r_m)³/12 = 125 005 mm⁴, is 3.2 % high: at t/r = 0.09 and with the web overlapping the arc's ends, the mid-line model is cruder than for the Z.

(c) Spar: 2 780 × 335.70 × 10⁻⁶ = **0.933 kg/m**. Rod: area π × 4.5² = 63.62 mm², 17 000 × 63.62 × 10⁻⁶ = **1.081 kg/m**. Combined centre of mass: (0.933 × (−12.646) + 1.081 × (−20.0))/(0.933 + 1.081) = **−16.59 mm**, that is **13.41 mm** behind the leading edge (x = −30 mm); the spar alone is at 17.35 mm.

(d) The centre of mass moves by (mass added) × (distance ahead of the current centre)/(total mass). A tungsten rod puts six times the mass of aluminium into a small volume, as far forward as it fits. Thickening the nose wall adds aluminium spread round the whole arc, whose centroid (−18.3 mm) lies behind the rod's position, so it would take more than six times the volume, and it changes the stiffness, which the rod barely does (it is not bonded to carry bending in this model).

**Check.** Symmetry about x gives ȳ = 0 and Ixy = 0 (the code returns 0.00). The full half-disc of radius 30 mm has its centroid at −4R/3π = −12.73 mm; the hollow spar's −12.65 mm is close, as expected for a section whose material lies near the same arc and diameter.

**Common mistakes.** Using the half-circle centroid 4R/3π for the thin arc (that is for the solid half-disc; the arc's is 2r/π); taking the inner arc down to x = 0, which deletes the web.

## 4. Your own Z-stringer (s = 15: h = 25 mm, b = 12 mm, t = 1.2 mm)

(a)

| method | A (mm²) | Ix (mm⁴) | Iy (mm⁴) | Ixy (mm⁴) |
|---|---|---|---|---|
| polygon (exact) | 55.92 | 5 236.1 | 1 188.7 | 1 850.69 |
| composite | 55.92 | 5 236.1 | 1 188.7 | 1 850.69 |
| thin-walled | 55.92 | 5 222.6 | 1 185.2 | 1 855.83 |

Thin-walled errors: Ix **−0.26 %**, Iy **−0.29 %**, Ixy **+0.28 %**.

(b) **I_max = 5 954.8 mm⁴, I_min = 470.0 mm⁴, θ_p = −21.22°**. `principal_axes_eig` returns eigenvalues 470.04 and 5 954.77; the eigenvector of 5 954.77 is (0.362, 0.932), at 68.78°: it is the *normal* of the I_max axis, so the axis is at 68.78° − 90° = −21.22°, as Mohr's circle says. `rotate_axes(Ix, Iy, Ixy, θ_p)` gives I_uv = 2 × 10⁻¹³ mm⁴.

(c) **r_min = √(470.0/55.92) = 2.90 mm**, against r_x = 9.68 mm: buckling as a column is governed by the inclined minor axis, not by either flange or web direction.

(d) The ratio I_min/I_max is largest at **b = 20.4 mm**, where it is 0.1155, with **θ_p = −40.2°** and **I_min = 1 481 mm⁴**. Wider flanges raise I_min, but beyond this the ratio falls because I_max grows faster.

**Check.** Ix + Iy = 6 424.8 = I_max + I_min; Ix Iy − Ixy² and I_max I_min are both 2.799 × 10⁶ mm⁸ (equal to round-off of the printed figures).

**Common mistakes.** Reporting the eigenvector's angle as the principal axis; using |Ixy| and so turning the axis the wrong way.

## 5. Idealise a wing box (s = 15: skins 1.5 mm)

(a) **Full Ix = 10 815 836 mm⁴**; booms only 6 048 000 mm⁴ (**−44.1 %**); lumped **10 944 000 mm⁴ (+1.18 %)**. Mass **9.32 kg/m** in 2024-T3.

(b) The two skins' own I_c + A d² about the centroid is 4 191 836 mm⁴: **38.8 %** of the full Ix. The skin is a major bending element, which is why dropping it is wrong.

(c)

| stringers per skin | full Ix (mm⁴) | lumped error | mass (kg/m) | Ix per kg/m (mm⁴ per kg/m) |
|---|---|---|---|---|
| 0 | 9 087 836 | +1.41 % | 7.99 | 1 137 441 |
| 2 | 9 951 836 | +1.29 % | 8.66 | 1 149 582 |
| 4 | 10 815 836 | +1.18 % | 9.32 | 1 159 985 |
| 8 | 12 543 836 | +1.02 % | 10.66 | 1 176 883 |

Each 60 mm² stringer pair adds 432 000 mm⁴ (2 × 60 × 60²). Ix per unit mass rises slowly, because stringers sit at the full skin height while the webs (which also carry mass) contribute little; the lumped error falls because the skin panels between booms get narrower, and the lumping rule is exact in the limit of many booms.

**Check.** Booms only does not depend on the skin thickness: 6 048 000 mm⁴ = 2 × (2 × 300 + 4 × 60) × 60² mm⁴.

**Common mistakes.** Lumping the skin into the booms *and* keeping it as rectangles (double counting).

## 6. How good are nominal section properties?

(a) The thin-walled Ix error of the 25/15/1.6 Z is −0.42 %. A thickness of 1.5 mm changes the exact Ix by **−5.15 %** and 1.7 mm by **+5.02 %** (area ±5.9 %). The tolerance matters about ten times more than the approximation; a margin calculation should use the minimum thickness, not the nominal.

(b) Two 4.8 mm holes through the 6 mm top flange (at x = ±15 mm) remove 57.6 mm². At that station the centroid moves **4.34 mm** down and **Ix falls by 8.4 %**; the top fibre, now farther from the centroid, sees more stress for two reasons. The spar's overall stiffness is barely changed, because the holes are short compared with the span; the net section matters for strength at the fastener row, the gross section for deflection.

(c) Fillets add material at the web–flange junction, near the flanges: they raise A and Ix slightly. For thin extrusions their effect is usually a fraction of a percent of Ix, smaller than the thickness tolerance.

(d) In preliminary design, nominal geometry with thin-walled formulas or lumped booms is adequate, with the minimum thickness if margins are tight. In a stress check at a fastener row, use the net section with minimum thickness; fillets can be ignored (conservatively) or included from the drawing.

## Check-yourself answers

1. x̄ = ∫x dA/A. For parts: x̄ = (Σ A_i x̄_i − A_h x̄_h)/(Σ A_i − A_h), with the hole's area and moment subtracted.
2. I = I_c + A d² and Ixy = Ixy_c + A d_x d_y; the cross term 2d ∫y_c dA vanishes only when y_c is measured from the centroid.
3. Second moment grows with distance squared; the flanges are far from the neutral axis (A d² large), while their own thickness is small (I_c ∝ t³, negligible: 900 mm⁴ against 1 555 200 mm⁴ in the lesson's spar).
4. A = ½ Σ (x_i y_{i+1} − x_{i+1} y_i). Listed clockwise, the same magnitude with a negative sign; all the integrals change sign together.
5. Each element x y dA has a mirror −x y dA (or x(−y) dA) across the axis of symmetry, so the integral cancels. A non-zero Ixy means x and y are not principal axes: bending about x also produces deflection along x.
6. I_max,min = (Ix + Iy)/2 ± √(((Ix − Iy)/2)² + Ixy²). Check: the eigenvalues of [[Iy, Ixy], [Ixy, Ix]] (`numpy.linalg.eigh`) are the same, and their eigenvectors are the normals of the principal axes.
7. About four times, 1.6 %: the error is of order (t/b)².
8. (947 × 2.114 + 20 × 3.00)/967 = **2.132 m**, that is 22.16 % MAC instead of 20.94 %.
