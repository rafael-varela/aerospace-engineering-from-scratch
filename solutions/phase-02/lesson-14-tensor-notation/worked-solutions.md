# P02.14  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. Your own fuselage seam ($s = 15$: p = 0.065 MPa, seam at 35°)

**Assumptions.** Thin-walled membrane, plane stress, axes 1 longitudinal, 2 hoop, 3 radial; MPa throughout.

**Derivation.** (a) σ₁₁ = pr/(2t) = 0.065 × 2.0/(2 × 0.0016) = 40.625 MPa, σ₂₂ = 81.25 MPa, σ₁₂ = 0.80 × 10⁶/(2π × 4.0 × 0.0016) = 19.894 MPa:

σ = [[40.625, 19.894, 0], [19.894, 81.25, 0], [0, 0, 0]] MPa.

(b) c = 60.9375, d = (40.625 − 81.25)/2 = −20.3125, cos 70° = 0.34202, sin 70° = 0.93969. σ_n = σ_y′ = 60.9375 + 20.3125 × 0.34202 − 19.894 × 0.93969 = **49.19 MPa**; τ_seam = |20.3125 × 0.93969 + 19.894 × 0.34202| = **25.89 MPa** (σ_x′ = 72.68 MPa along the seam). `seam_stresses`, the rotated tensor's (2,2) and (1,2) entries and `mohr_2d` all give 49.190 and 25.892 MPa.

(c) I₁ = 121.875 MPa, I₂ = 40.625 × 81.25 − 19.894² = 2904.99 MPa², I₃ = 0. R = √(20.3125² + 19.894²) = 28.432 MPa; σ₁ = 60.9375 + 28.432 = **89.37 MPa**, σ₂ = **32.51 MPa**, σ₃ = 0. tan 2θ_p = 2 × 19.894/(40.625 − 81.25) = −0.9794 with σ_x − σ_y < 0, so 2θ_p = 135.60°, θ_p = **67.80°** from the longitudinal axis. `principal_stresses`: 89.370, 32.505, 0 MPa and the first column at 67.80°.

(d) τ_max = σ₁/2 = **44.68 MPa** (not R = 28.43: the plane is through the thickness). σ_vm = √(89.37² − 89.37 × 32.51 + 32.51²) = **78.35 MPa**. Margin = 325/78.35 − 1 = **3.15**.

**Check.** σ₁ + σ₂ = 121.88 = I₁ ✓; σ₁σ₂ = 2905.0 = I₂ ✓. The seam normal stress 49.19 lies between σ₂ and σ₁ ✓, and at a = 0 the code gives 81.25 (hoop) ✓.

**Common mistakes.** Using σ_x′ (72.68) for the seam instead of σ_y′; quoting R as the maximum shear; taking θ_p = ½ atan(−0.9794) = −22.2°, which is the σ₂ direction (67.80° − 90°).

## 2. Index gymnastics

(a) Fix i. ε_ijk ≠ 0 needs {j, k} = the two indices other than i, in some order; ε_imn likewise for {m, n}. If {j, k} ≠ {m, n} some δ on the right vanishes pairwise and the left is zero. If (j, k) = (m, n) the left is ε_ijk² = 1 and the right is 1 − 0. If (j, k) = (n, m) the left is −1 and the right is 0 − 1. Nine cases, all agree. Then (a × b)·(c × d) = ε_ijk a_j b_k ε_imn c_m d_n = (δ_jm δ_kn − δ_jn δ_km) a_j b_k c_m d_n = (a·c)(b·d) − (a·d)(b·c). With c = a, d = b: |a × b|² = |a|²|b|² − (a·b)².

(b) A passive rotation of the axes by a about n carries e_i to e′_i = cos a e_i + (1 − cos a)(e_i · n) n + sin a (n × e_i) (Rodrigues for a vector: keep the part along n, rotate the perpendicular part by a in the plane spanned by e_i⊥ and n × e_i). Take the component along e_j: Q_ij = e′_i · e_j = cos a δ_ij + (1 − cos a) n_i n_j + sin a (n × e_i)_j. The last: (n × e_i)_j = ε_jkl n_k (e_i)_l = ε_jki n_k = ε_ijk n_k by cyclic permutation. This is (3). Check: for n = e₃, a = 90°: Q = [[0, 1, 0], [−1, 0, 0], [0, 0, 1]], so e₁ has new components (0, −1, 0) ✓.

(c) In the old frame t_p = σ_pq n_q. Both t and n are vectors: t′_i = Q_ip t_p and n′_j = Q_jq n_q. Invert the second: multiply by Q_jr and sum on j, Q_jr n′_j = Q_jr Q_jq n_q = δ_rq n_q = n_r, **using orthogonality QᵀQ = I**. Then t′_i = Q_ip σ_pq n_q = Q_ip σ_pq Q_jq n′_j = (Q_ip Q_jq σ_pq) n′_j. Cauchy's formula holds in the new frame, t′_i = σ′_ij n′_j, for every n′, so σ′_ij = Q_ip Q_jq σ_pq.

(d) σ′_ii = Q_ip Q_iq σ_pq = δ_pq σ_pq = σ_pp. σ′_ij σ′_ji = Q_ip Q_jq σ_pq Q_jr Q_is σ_rs = δ_ps δ_qr σ_pq σ_rs = σ_pq σ_qp. det σ′ = det Q · det σ · det Qᵀ = det σ. Rotation about e₃: Q = [[c, s, 0], [−s, c, 0], [0, 0, 1]]. σ′₁₁ = Q₁p Q₁q σ_pq = c²σ_x + cs τ + sc τ + s²σ_y = c²σ_x + s²σ_y + 2cs τ; σ′₂₂ = s²σ_x + c²σ_y − 2cs τ; σ′₁₂ = Q₁p Q₂q σ_pq = −cs σ_x + c² τ − s² τ + sc σ_y = −(σ_x − σ_y) cs + (c² − s²) τ. With c² = (1 + cos 2θ)/2, s² = (1 − cos 2θ)/2, 2cs = sin 2θ, c² − s² = cos 2θ these are (9).

**Common mistakes.** Writing ε_jik for the cross term and getting the active rotation; forgetting that inverting (4) needs orthogonality, not just Q⁻¹ (the law (5) is Q σ Qᵀ, not Q σ Q⁻¹ in general, and the two coincide only because Q is orthogonal).

## 3. A strain-gauge rosette ($s = 15$: ε_a = 950, ε_b = 300, ε_c = −200 microstrain)

**Assumptions.** Small strain, isotropic linear elastic, plane stress, gauges exactly at 0°, 45°, 90°; microstrain (10⁻⁶) throughout.

**Derivation.** (a) ε_x = ε_a = 950, ε_y = ε_c = −200. At θ = 45°, (9) for strain gives ε_x′ = (ε_x + ε_y)/2 + ε_xy, so ε_xy = ε_b − (ε_a + ε_c)/2 = 300 − 375 = **−75 microstrain** (γ_xy = −150).

(b) Centre 375, R = √(575² + 75²) = 579.87: ε₁ = **954.9**, ε₂ = **−204.9 microstrain**. tan 2θ_p = 2(−75)/1150 = −0.1304, 2θ_p = −7.43°, θ_p = **−3.72°** (the code reports 176.28°, the same line). `principal_stresses` on the strain tensor: 954.87, −204.87 and, once ε₃₃ is included, −369.40.

(c) Inverse law with σ₃₃ = 0: ε₃₃ = −ν(σ₁₁ + σ₂₂)/E. With the plane-stress formulas σ₁₁ + σ₂₂ = E(ε_x + ε_y)/(1 − ν), so ε₃₃ = −ν(ε_x + ε_y)/(1 − ν) = −0.33 × 750/0.67 = **−369.4 microstrain**. The full tensor is ε = [[950, −75, 0], [−75, −200, 0], [0, 0, −369.4]] × 10⁻⁶. `hooke_stress(isotropic_stiffness(73e9, 0.33), eps)` gives

σ = [[72.42, −4.12, 0], [−4.12, 9.30, 0], [0, 0, 0.00]] MPa,

with σ₃₃ = −7 × 10⁻⁹ Pa. Plane-stress formulas: σ_x = 73 × 10³ (950 − 0.33 × 200) × 10⁻⁶/(1 − 0.1089) = 72.42 MPa ✓, σ_y = 73 × 10³(−200 + 0.33 × 950) × 10⁻⁶/0.8911 = 9.30 MPa ✓, τ_xy = 73 × 10³ × (−75) × 10⁻⁶/1.33 = −4.12 MPa ✓. Principal stresses 72.69 and 9.03 MPa on the same axes as the principal strains (isotropy); σ_vm = 68.6 MPa, margin 3.74.

(d) With ε₃₃ = 0 the code gives [[92.10, −4.12, 0], [−4.12, 28.98, 0], [0, 0, 39.95]] MPa: every diagonal term gains λ × 369.4 × 10⁻⁶ = 19.68 MPa (λ = 52.3 GPa), and a **spurious σ₃₃ = 39.95 MPa** appears on the free surface, where the stress must be zero. The shear is unaffected.

**Check.** ε₁ + ε₂ = 750 = ε_x + ε_y ✓. The strains are consistent with a skin stressed mainly along the gauge-a direction, and the free-surface condition σ₃₃ = 0 is recovered to 10⁻⁹ Pa by the 3-D tensor law once ε₃₃ is correct ✓.

**Common mistakes.** Using γ_xy = 2ε_xy in the tensor; forgetting ε₃₃ (part (d) shows the size of the error: 27 % on σ_x); assuming plane strain (ε₃₃ = 0) for a thin sheet, which is plane stress.

## 4. Seam angle and torque ($s = 15$)

**Derivation.** (a) Sweep at 0.1°: maximum normal stress **89.37 MPa at 157.8°**, minimum **32.51 MPa at 67.8°**, maximum shear **28.43 MPa at 22.8°** (and 112.8°). The seam with the largest opening stress is perpendicular to σ₁ (θ_p + 90° = 157.8°), the seam with the smallest is along σ₁'s direction (67.8°, its normal is the σ₂ direction), both carry zero shear, and the shear maximum R = 28.43 MPa sits 45° from them. Your 35° seam is 12° from the worst shear angle.

(b) With T → −T, the 35° seam sees **86.58 MPa normal and 12.28 MPa shear**: the opening stress rises from 49.2 to 86.6 MPa. Under the original torque the seam at −35° (= 145°) sees the same. The invariants and σ_vm (78.35 MPa) are unchanged.

(c) `brentq` on σ_vm − 325/1.5 MPa over (0, 50) MN m: **T* = 4.758 MN m**, where σ₁₂ = 118.3 MPa. By hand: σ_vm² = σ_x² − σ_xσ_y + σ_y² + 3τ²; (216.67)² = 46944; σ_x² − σ_xσ_y + σ_y² = 1650.4 − 3300.8 + 6601.6 = 4951.2; τ* = √((46944 − 4951)/3) = 118.3 MPa; T* = 2π × 4.0 × 0.0016 × 118.3 × 10⁶ = 4.758 MN m ✓. Seam stresses at T*: −43.3 MPa (the seam is *closed*) and 59.6 MPa shear.

(d) From the code (MPa, degrees):

| T [MN m] | σ₁₂ | seam normal | seam shear | σ₁ | σ₂ | σ₃ | θ_p | σ_vm | margin |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0 | 67.88 | 19.09 | 81.25 | 40.62 | 0 | 90.0 | 70.36 | 3.62 |
| 0.8 | 19.89 | 49.19 | 25.89 | 89.37 | 32.51 | 0 | 67.8 | 78.35 | 3.15 |
| 1.6 | 39.79 | 30.50 | 32.70 | 105.61 | 16.26 | 0 | 58.5 | 98.49 | 2.30 |
| 2.4 | 59.68 | 11.80 | 39.50 | 123.98 | 0 | −2.1 | 54.4 | 125.05 | 1.60 |
| 3.2 | 79.58 | −6.89 | 46.30 | 143.07 | 0 | −21.2 | 52.2 | 154.75 | 1.10 |
| 4.0 | 99.47 | −25.59 | 53.11 | 162.46 | 0 | −40.6 | 50.8 | 186.11 | 0.75 |

The in-plane second principal stress crosses zero when I₂ = σ_xσ_y − τ² = 0, i.e. τ = √(40.625 × 81.25) = 57.5 MPa, T = 2.31 MN m; above that the skin carries a diagonal compression, and a thin panel between frames and stringers fails by **shear buckling** (a diagonal-tension field) long before yield. The yield margin, which is a function of invariants only, is blind to this.

**Check.** At T = 0 the hoop stress is principal (θ_p = 90°) ✓ and the seam normal stress 67.88 = c − d cos 70° = 60.94 + 20.31 × 0.342 ✓. The T* found numerically equals the hand value to four figures ✓.

**Common mistakes.** Reading the seam angle of maximum normal stress as θ_p instead of θ_p + 90°; stopping at two principal stresses and missing the negative third one; treating the margin as a complete safety statement for a thin panel.

## 5. A composite lamina as a fourth-order tensor ($s = 15$: θ = 45°)

**Assumptions.** Plane stress of a single lamina; the fourth-order tensor holds only the in-plane entries; GPa throughout.

**Derivation.** ν₂₁ = 0.30 × 10/140 = 0.021429; 1 − ν₁₂ν₂₁ = 0.993571. Q₁₁ = **140.906**, Q₂₂ = **10.065**, Q₁₂ = **3.019**, Q₆₆ = **5.000 GPa**.

(a) Filling all four index orders of (12) gives a tensor that passes `np.allclose` for C_ijkl = C_jikl, C_ijlk and C_klij ✓. `hooke_stress(C, diag(ε, 0, 0))` returns σ₁₁ = Q₁₁ε, σ₂₂ = Q₁₂ε, σ₁₂ = 0 ✓.

(b) `transform_tensor4(C, rotation_matrix([0, 0, 1], -45°))`: C̄₁₁₁₁ = **44.252**, C̄₂₂₂₂ = **44.252**, C̄₁₁₂₂ = **34.252**, C̄₁₂₁₂ = **36.233**, C̄₁₁₁₂ = **32.710 GPa** (and C̄₂₂₁₂ = 32.710). Classical formulas with c = s = 1/√2 (c⁴ = s⁴ = c²s² = ¼): Q̄₁₁ = ¼(140.906 + 2(3.019 + 10) + 10.065) = 44.252 ✓; Q̄₁₆ = ¼(140.906 − 3.019 − 10) − ¼(10.065 − 3.019 − 10) = 31.972 + 0.739 = 32.710 ✓. At 45° the lamina is symmetric in x and y, which is why C̄₁₁₁₁ = C̄₂₂₂₂ and C̄₁₁₁₂ = C̄₂₂₁₂. (Rotating the axes by +45° instead gives C̄₁₁₁₂ = −32.710: the fibres would then be at −45°.)

(c) Q̄ = [[44.252, 34.252, 65.421], [34.252, 44.252, 65.421], [32.710, 32.710, 72.466]] GPa (third column doubled). Inverting: **E_x = 13.21 GPa, E_y = 13.21 GPa, G_xy = 17.95 GPa, ν_xy = 0.321**, and η_xs = S̄₁₃/S̄₁₁ = −0.613. For σ_x = 100 MPa: ε_x = **7571**, ε_y = **−2429**, γ_xy = **−2321 microstrain**. The shear strain is non-zero because the compliance S̄₁₃ ≠ 0: an off-axis lamina couples extension to shear (the fibres try to straighten towards the load), and only at 0° and 90° does the coupling vanish. Note G_xy = 17.9 GPa at 45°, far above G₁₂ = 5 GPa: a ±45° layer is how shear-stiff skins are made.

(d) E_x [GPa]: 0°: 140.0; 15°: 52.9; 30°: 21.3; 45°: 13.2; 60°: 10.7; 75°: 10.1; 90°: 10.0. E_x = 14 GPa (one-tenth of E₁) between 30° and 45°, at about **43°** by interpolation on 1/E_x. Half the stiffness is gone by 15°: fibre alignment tolerance matters.

**Check.** At θ = 0 the inverse of Q̄ returns E_x = 140.0, E_y = 10.0, G_xy = 5.0, ν_xy = 0.300 exactly ✓. The table is symmetric under θ → 90° − θ with E_x ↔ E_y ✓.

**Common mistakes.** Filling only C₁₂₁₂ (the tensor then lacks the minor symmetry and the rotation is wrong); forgetting the factor 2 on the third column (γ = 2ε); the sign of θ, which flips C̄₁₁₁₂ and η_xs but nothing else.

## 6. Is the seam acceptable? ($s = 15$)

(a) Von Mises: 325/78.35 − 1 = **3.15**. Tresca: 325/(2 × 44.68) − 1 = **2.64**. The ratio σ_vm/(2τ_max) is 0.877 here, within the allowed band [0.866, 1]. For a riveted seam, whose rivets are loaded in shear and whose sheet fails by bearing and shear-out at the holes, quote the **Tresca** (maximum-shear) margin: it is never less conservative, and the difference (16 %) is of the order of the uncertainty in the load.

(b) Ranked by effect on the local normal stress at the seam: (1) stress concentration at the rivet holes, factor ≈ 3 upward at the hole edge, 49 → 148 MPa; (2) secondary bending at the lap, a further factor ≈ 1.5 on the outer surface of one sheet, 148 → 220 MPa, approaching yield locally; (3) frame and stringer load sharing, downward by tens of percent on the longitudinal and shear components (which lowers the seam stresses by less, since the hoop stress dominates them); (4) the thin-wall assumption, exact to t/r = 0.0008, negligible.

(c) The seam is loaded by the pressurisation cycle: every flight applies and removes 49 MPa across it on average and about 150–220 MPa at the hole edge. That stress *range*, repeated tens of thousands of times, is what sizes a fuselage seam, through the fatigue curve of the riveted joint, not the static yield margin. A margin of 3 on yield says nothing about crack initiation at 20 000 cycles, which is why fuselage seams are the classic site of fatigue failures and inspections.

(d) The principal *stresses* 89.37, 32.51 and 0 MPa are invariants: the same in every frame, as are I₁, I₂, I₃, τ_max, σ_vm and both margins. What depends on the frame is the list of components σ_ij and the *components* of the principal directions; the directions themselves are physical lines in the skin, and the angle 67.8° is measured from the fuselage axis only because we chose to quote it that way. The seam stresses are neither: they belong to a specific plane, and the frame only chooses how we compute them. So compute in whatever frame is easiest and quote the invariants.

**Check.** The Tresca margin equals the von Mises margin only for uniaxial stress; here σ₂/σ₁ = 0.36, so they must differ, and by less than 15.5 % ✓.

**Common mistakes.** Quoting the in-plane R for Tresca (gives a non-conservative 325/56.9 − 1 = 4.7); treating a yield margin as a fatigue margin.
