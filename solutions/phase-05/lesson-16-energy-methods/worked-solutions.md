# P05.16  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. Dummy loads on a cantilever (s = 15: q = 2.5 kN/m)

**Assumptions.** Linear elastic, slender, bending energy only; x from the root.

(a) M = q(L − x)²/2 + Q(L − x). ∂U/∂Q|_{Q=0} = ∫₀ᴸ [q(L − x)²/2](L − x)/EI dx = qL⁴/(8EI) = 2500 × 16/(8 × 50 000) = **0.100 m = 100 mm**.

(b) M = q(L − x)²/2 + C. ∂U/∂C|_{C=0} = ∫₀ᴸ q(L − x)²/(2EI) dx = qL³/(6EI) = **0.0667 rad (66.7 mrad)**.

(c) m = L/2 − x on 0 ≤ x ≤ L/2. δ = ∫₀^{L/2} q(L − x)²(L/2 − x)/(2EI) dx = **35.42 mm**, equal to q x²(6L² − 4Lx + x²)/(24EI) at x = 1 m: 2500 × 1 × (24 − 8 + 1)/(24 × 50 000) = 35.42 mm.

(d) `castigliano_deflection` on U(Q, C) gives **100.000 mm** and **66.667 mrad**, to 10⁻⁹ relative.

**Check.** Units: N/m × m⁴/(N·m²) = m. The tip slope times 3L/4 equals the tip deflection for a uniform load (66.7 mrad × 1.5 m = 100 mm), a standard result.

**Common mistakes.** Setting Q = 0 before differentiating; measuring x from the tip in one term and from the root in another.

## 2. Derivations

(a) With δ = FP, Clapeyron gives U = ½ΣP_kδ_k = ½ΣΣF_jkP_jP_k. Then ∂U/∂P_j = ½Σ_k(F_jk + F_kj)P_k = Σ_kF_jkP_k = δ_j, using F_jk = F_kj from (b). **δ_j = ∂U/∂P_j.**

(b) Apply P_a first: work ½F_aaP_a². Then P_b: work ½F_bbP_b² by P_b, plus F_abP_b × P_a by P_a, which is already at full value while its point moves F_abP_b. Total ½F_aaP_a² + ½F_bbP_b² + F_abP_aP_b. In the other order the cross term is F_baP_bP_a. The final state is the same and an elastic body's energy depends only on its state, so **F_ab = F_ba**: the deflection at a due to a unit load at b equals that at b due to a unit load at a.

(c) Release the prop: M = q(L − x)²/2 − R(L − x). ∂U/∂R = −∫₀ᴸ M(L − x)/EI dx = −(qL⁴/8 − RL³/3)/EI = 0 ⇒ **R = 3qL/8**. Mid-span deflection: the real M with R = 3qL/8 and, as the virtual system, a unit load at L/2 on the released cantilever (it is in equilibrium, which is all the principle needs): δ = ∫₀^{L/2} M (L/2 − x)/EI dx = **qL⁴/(192EI)** (checked symbolically). The closed form of the propped beam, q x²(3L² − 5Lx + 2x²)/(48EI) at L/2, gives the same.

(d) As in the lesson: M = M_q − S sin θ (a − x) on x < a, N_spar = −S cos θ on x < a, N_strut = S. U = ∫M²/2EI + (S cos θ)²a/(2EA_spar) + S²L_s/(2EA_strut). ∂U/∂S = −sin θ∫₀ᵃ M(a − x)/EI dx + S cos²θ a/EA_spar + S L_s/EA_strut = 0. Substitute M and collect S: **S = sin θ ∫₀ᵃ M_q(a − x)/EI dx / (sin²θ ∫₀ᵃ(a − x)²/EI dx + cos²θ a/EA_spar + L_s/EA_strut).**

**Common mistakes.** In (b), giving the cross term a factor ½ (the first load does not rise during the second loading); in (c), using a virtual moment that is not in equilibrium with the unit load.

## 3. A curved frame segment (s = 15: R = 0.55 m)

**Assumptions.** Thin curved bar (radius large compared with depth), bending energy only.

(a) At angle φ the point is (R cos φ, R sin φ); the vertical load at B = (0, R) has lever arm R cos φ: **M = PR cos φ**. U = ∫₀^{π/2} P²R²cos²φ R dφ/(2EI) = **πP²R³/(8EI)**. **δ_v = πPR³/(4EI) = π × 2000 × 0.166375/(4 × 15 000) = 17.42 mm.**

(b) A horizontal unit load at B has lever arm equal to the vertical distance R − R sin φ: m = R(1 − sin φ). δ_h = ∫ PR cos φ · R(1 − sin φ) R dφ/EI = (PR³/EI)(1 − ½) = **PR³/(2EI) = 11.09 mm**; with the signs worked through (the free part of the arc rotates in the sense of the load's moment), B moves away from A, in the −x direction, as well as down.

(c) A straight cantilever of length R: PR³/(3EI) = **7.39 mm**. The curved bracket is 2.36 times more flexible: its arc is π/2 times longer and the moment stays near its maximum PR along much of the arc near A.

(d) δ_shear/δ_bend = (αPL/GA)/(PL³/3EI) = **3αEI/(GAL²)**. For 25 × 150 mm: I/A = h²/12 = 1.875 × 10⁻³ m². Ratio = 3 × 1.2 × 73.1 × 1.875 × 10⁻³/(28 L²): **0.196 at 0.3 m, 0.0176 at 1.0 m, 0.0007 at 5.0 m**. Neglecting shear is unconservative by about 20 % for a short, deep beam (L/h = 2) and negligible for a slender spar (L/h > 10); it matters most for thin-webbed sections, where A is the web area only.

**Check.** Numerical Simpson integration of the two integrals gives 17.4227 and 11.0917 mm.

**Common mistakes.** Forgetting the R in the arc length ds = R dφ; using the straight distance instead of the lever arm.

## 4. The engine-mount truss under your own load (s = 15: 6.0 kN at −60°)

**Assumptions.** Pin joints, small deflections, all three bars the same 4130 tube (A = 50.72 mm², EA = 10.40 MN).

(a) F = (3.000, −5.196) kN. Bar forces (bottom, middle, top), tension positive: **−3.890, +1.807, +5.714 kN**, identical from the force and stiffness methods.

(b) Unit load along F: **δ = 0.4813 mm**; the stiffness displacement (0.0342, −0.5360) mm projected on the load direction: **0.4813 mm**.

(c) Reactions (from (C N) at the support degrees of freedom): bottom wall node **(3.479, 1.740) kN**, middle **(−1.798, 0.180) kN**, top **(−4.681, 3.277) kN**. Sum (−3.000, +5.196) kN: equal and opposite to F, ΣF = 0 to round-off. Moments about the origin: the wall nodes have x = 0, so their reactions contribute −yR_x = 0.540 + 2.809 = 3.348 kN·m; the load contributes 0.5 × (−5.196) − 0.25 × 3.000 = −3.348 kN·m. **ΣM = 0.**

(d)

| middle-bar EA factor | bottom [kN] | middle [kN] | top [kN] |
|---|---|---|---|
| 0.1 | −3.052 | +0.300 | +6.629 |
| 1 | −3.890 | +1.807 | +5.714 |
| 10 | −4.904 | +3.631 | +4.606 |
| 100 | −5.131 | +4.038 | +4.359 |

The redundant bar's share rises with its stiffness and saturates: as its EA → ∞ it becomes inextensible, the engine node can then move only perpendicular to it, and the forces tend to a limit set by how the other two bars resist that motion. Stiffness attracts load, with diminishing returns.

**Check.** At factor 1 the unit-load displacement equals the stiffness displacement; ΣF and ΣM vanish.

**Common mistakes.** Building C on the free degrees of freedom only and then trying to read reactions from it; mixing the sign convention (tension positive) between methods.

## 5. Where should the strut go? (s = 15: h = 1.10 m)

**Assumptions.** As the lesson; 1001 stations in the sweep.

(a) Elliptic lift: the largest bending moment is smallest at **a = 3.00 m**, max |M| = 4.58 kN·m, S = 33.78 kN, R = 11.63 kN, tip deflection 53.1 mm.

(b) Uniform lift: best **a = 3.40 m**, max |M| = 4.93 kN·m, S = 39.72 kN, tip deflection 59.5 mm. The uniform distribution puts more lift outboard (its centroid is at L/2 = 2.5 m against 4L/(3π) = 2.12 m for the ellipse), so the overhang carries more moment and the strut must move outboard to balance the root and strut-point moments.

(c) At a = 3.0 m: elliptic S = 33.78 kN (R = 11.63 kN), root moment 4.58 kN·m; uniform S = 42.52 kN (R = 14.64 kN), root moment 2.58 kN·m, tip 92.7 mm. The lower strut (h = 1.1 m instead of 1.2 m) makes the strut flatter, so the same vertical reaction needs more tension: compare the lesson's 31.5 kN at h = 1.2 m.

**Check.** At each optimum the root moment and the moment at the strut point are close (4.58 against 4.32 kN·m elliptic, 4.93 against 4.76 kN·m uniform, on a 0.05 m grid): moving the strut either way would raise one of them, which is the condition for the smallest maximum.

**Common mistakes.** Minimising the root moment instead of the largest moment anywhere (the root moment can be driven negative, overloading the strut station).

## 6. Judgment: the finite-element model disagrees

A model plan, in outline. **Candidate differences and their direction.** (1) Strut-end fittings or bolts modelled with flexibility: a softer strut path, S falls. (2) Spar root modelled as a flexible fitting rather than built in: the spar root rotates, the spar's deflection at the strut point grows, S rises; a pinned root would raise it a lot. (3) Spar shear deformation included: the spar is more flexible, S rises slightly. (4) Strut attachment offset from the spar's neutral axis, or a different h: sin θ changes; a flatter strut needs more tension for the same vertical component, but carries a smaller share. (5) Lift applied at the quarter-chord with the wing's torsion modelled, or a different spanwise distribution (more inboard lift lowers S). (6) Geometric non-linearity switched on. (7) A units error in EA or a wrong cross-section. **First checks.** The load-path totals: the FE strut's vertical component plus the FE root shear must equal 18.6 kN, and the root moment must follow from the lift moment minus R × a. Then the FE strut's axial stiffness EA/L_s against the hand value. **Settling measurement.** Strain gauges on the strut (and on the spar caps at the root) during a ground proof-load test with distributed sandbags: measured strut tension against both predictions.

## Check-yourself answers

1. **U = N²L/(2EA)**, **U = T²L/(2GJ)**, **U = ∫M²/(2EI) dx.**
2. **δ_j = ∂U/∂P_j**, the displacement at a load in its direction. Add a dummy load Q at the point and in the direction wanted, differentiate U with respect to Q, then set Q = 0.
3. **δ = ∫ M m/EI dx**, where m(x) is the bending moment produced by a unit load applied at the point and in the direction of the wanted deflection.
4. The deflection at point a due to a unit load at b equals the deflection at b due to a unit load at a; the **flexibility matrix is symmetric**.
5. **R = 3qL/8**, from **∂U/∂R = 0** (the prop does not move: compatibility).
6. The tip load also **twists** the leg with torque Pa, adding Pa²b/(GJ), which here is more than half of the total.
7. **S → 0**: the soft strut stretches freely, so it cannot resist the spar's deflection; in (8) the term L_s/EA_strut dominates the denominator.
8. The **redundant forces** (one unknown per degree of indeterminacy). It is the smaller problem when the structure has few redundants and many joints; the stiffness method needs one unknown per free degree of freedom.
