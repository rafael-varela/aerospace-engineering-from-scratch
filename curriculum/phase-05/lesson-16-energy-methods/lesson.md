# P05.16  Energy Methods

> A loaded structure stores the work done on it. Ask how fast that energy grows with a load and you have the deflection under it; ask which redundant force makes it stationary and you have solved an indeterminate structure.

| | |
|---|---|
| **Learning outcomes** | 1. Write the strain energy of members in tension, bending and torsion, and use Castigliano's second theorem δ_j = ∂U/∂P_j to find a deflection, introducing a dummy load where none acts. 2. Express the same result as the unit-load (virtual-force) method, assemble it as a flexibility matrix, and show that the matrix is symmetric (Maxwell's reciprocal theorem) and that its quadratic form is the energy (Clapeyron). 3. Solve a statically indeterminate structure by the force method: choose the redundant, impose ∂U/∂X = 0, and check the result by compatibility and by a limiting case. 4. Implement and test this on a spar under elliptic lift, an L-bracket in bending and torsion, an engine-mount truss and a strut-braced wing, with direct integration and the direct stiffness method as independent checks. |
| **Prerequisites** | `P05.03` (trusses, method of joints), `P05.07` (axial deformation NL/EA), `P05.08` (torsion, TL/GJ), `P05.09` and `P05.10` (bending moment diagrams, flexure), `P05.12` (deflection by integration and superposition), `P01.10` (integrals), `P01.13` (linear systems and null spaces), `P03.08` (direct solvers), `P03.04` (quadrature: Simpson's rule). Skill check: you can write the bending moment along a cantilever under a tip load and under a uniform load. |
| **Engineering question** | A light aircraft's wing has a spar 5 m long from the fuselage to the tip (2024-T3, EI = 292 kN·m²) carrying 18.6 kN of elliptically distributed lift at limit load, and a lift strut from the bottom of the fuselage, 1.2 m below the spar, to the spar 3 m out. The strut makes the wing statically indeterminate. How much load does the strut carry, and what does it do to the root bending moment and the tip deflection? The answers: the strut is in tension at 31.5 kN, lifting 11.7 kN of the 18.6 kN off the spar; the root moment falls from 39.5 kN·m (unbraced) to 4.4 kN·m and the tip deflection from 769 mm (an unbraced spar this light would be outside small-deflection theory) to 49 mm. |
| **Assumptions and limits** | Linear elastic members and small deflections, so the energy is a quadratic function of the loads and superposition holds; slender members, so the shear-deformation energy in bending is neglected (it is a few per cent for a deep spar and is estimated in Problem 3); pin-jointed trusses; for the braced wing, a spar built in at the fuselage, a strut pinned at both ends, and no change of geometry under load (the strut angle stays fixed). Energy methods give displacements at chosen points and redundant forces; they do not give stresses directly, and they say nothing about stability (`P05.15`). |
| **Analytical work** | Strain energy N²L/2EA, ∫M²/2EI dx, T²L/2GJ; Castigliano's second theorem and the dummy-load trick; the unit-load integral ∫Mm/EI dx; the tip deflection of a cantilever under a tip load, a uniform load and elliptic lift; the L-bracket's bending-plus-torsion deflection; Maxwell's reciprocal theorem; the propped cantilever's reaction 3qL/8 from ∂U/∂R = 0; the braced spar's redundant strut force with spar bending, spar axial and strut axial energy. |
| **Computational deliverable** | `starter/main.py`: `strain_energy_axial`, `strain_energy_torsion`, `strain_energy_bending`, `castigliano_deflection`, `unit_load_deflection`, `flexibility_matrix`, `strain_energy_members`, `l_bracket_members`, `l_bracket_tip_deflection`, `cantilever_moment`, `deflection_by_integration`, `truss_equilibrium_matrix`, `truss_force_method`, `truss_unit_load_displacement`, `truss_stiffness_method`, `braced_spar` (the `Member` dataclass, `elliptic_lift`, `truss_lengths` and `main` are given). Tests prove: the energy formulas; Castigliano's numerical ∂U/∂P equals PL³/3EI to 10⁻⁹; the unit-load tip deflection equals direct double integration of M/EI for three lift distributions; the L-bracket matches its closed form and its torsion-free limit; the flexibility matrix is symmetric, positive definite and gives U = ½PᵀFP; superposition; the force method equals the direct stiffness method on a determinate and an indeterminate truss, with joint equilibrium to 10⁻⁸; a vanishing redundant bar carries nothing; the braced spar satisfies compatibility and global equilibrium, recovers 3qL/8 when rigid, and its strut force grows monotonically with strut stiffness from zero. |
| **Practical activity** | **computational experiment**: unit-load against direct integration for three lift distributions; Castigliano on the L-bracket term by term; force method against stiffness method on the engine-mount truss with and without a redundant bar; a sweep of strut stiffness from zero to rigid and of strut position along the span. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `energy_methods.png`, the hand derivations of the L-bracket deflection and the propped-cantilever reaction, the experiment tables, and a half-page note on where to put the lift strut. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: one deflection and one redundant are derived by hand with the energy method; all public tests pass; the braced-spar result is checked by compatibility and by a limiting case; the note justifies the strut position from the sweep. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Roylance, modules "Trusses": Energy approach; "Beam Displacements": Energy method; "Shear and Torsion": Energy method for rotational displacement (required-free): strain energy and Castigliano on the three member types. `@libretexts-wierzbicki-structural-mechanics` Ch. 7 "Energy Methods in Elasticity": 7.1 The Concept of Potential Energy; 7.2 Equivalence of the Minimum Potential Energy and Principle of Virtual Work; 7.6 Castigliano Theorem (supplementary-free): the variational view behind the theorems. |
| **Extensions** | ★ Add the shear-deformation energy ∫αV²/(2GA) dx to the spar and show when it matters (deep, short beams). ★ Make the strut angle change with the deflection and iterate; how large must the lift be before the geometric change matters? ★ Write the braced spar with two struts (a V-strut) as a two-redundant problem: the force method then needs a 2 × 2 compatibility system; check it against a finite-difference beam with two spring supports. |
| **Time** | unmeasured |

## 1. The engineering problem

Many light aircraft carry their wing on a lift strut: a streamlined tube from the bottom of the fuselage to the wing spar a little more than halfway out. Without it the spar is a cantilever, and its root bending moment, the number that sizes the spar caps and the wing-attach fittings, is the full moment of the lift about the root. With it, part of the lift goes down the strut into the fuselage, the spar becomes a beam on two supports with an overhang, and its bending moment drops by an order of magnitude. That is why a strut-braced wing can be so light. But how much load goes into the strut? Equilibrium alone cannot say: the root of a cantilever spar can already carry the whole lift, so the strut is a **redundant** member, and its force depends on how stiff it is compared with the spar. Find it wrong and you either size the strut for too little load or the spar for too little moment.

The tools for this are the **energy methods**. An elastic structure under load stores strain energy U equal to the work done on it. Castigliano's second theorem says that the deflection under a load, in its direction, is ∂U/∂P. The same idea in reverse says that a redundant force X must take the value that makes U stationary, ∂U/∂X = 0, which is the statement that the structure fits together (compatibility). Written as the **unit-load method**, the deflection is the integral of the real internal forces times those of a unit load; written as a matrix over all loads, it is the **flexibility matrix**, whose symmetry is Maxwell's reciprocal theorem. One framework handles axial, bending and torsion members, trusses and beams, determinate and indeterminate structures alike.

Success criteria, written first: the unit-load deflection equals direct integration of M/EI to the quadrature error; Castigliano's numerical derivative equals the closed forms; the force method equals the direct stiffness method on a truss; the braced spar satisfies compatibility at the strut to 10⁻⁴ and gives 3qL/8 in its rigid limit.

## 2. Understand and derive

### 2.1 Strain energy in the three member types

A bar of length L carrying a constant axial force N stretches by NL/EA; the force rises from 0 to N as the stretch grows, so the work done, stored as energy, is half force times stretch:

**U_axial = N²L / (2EA).**   (1)

A shaft carrying torque T twists by TL/GJ (`P05.08`), and in the same way **U_torsion = T²L/(2GJ)**. A beam element dx under bending moment M rotates by dθ = M dx/EI, so

**U_bending = ∫ M² / (2EI) dx.**   (2)

The energies of the members of a structure add. Shear force V also stores energy, ∫αV²/(2GA) dx with α a shape factor; for slender beams it is small compared with (2) and is neglected here.

### 2.2 Castigliano's second theorem

Load a linear elastic structure with forces P₁ … P_n and let δ_j be the displacement at P_j, in its direction. The work done is U = ½ΣP_jδ_j (**Clapeyron's theorem**: each load rises linearly with its displacement). Now increase P_j alone by dP_j. The extra energy stored is dU; if instead you had applied dP_j first and then all the loads, the dP_j would have been carried through the full displacement δ_j, doing work dP_j δ_j. Linearity makes the final state the same either way, so

**δ_j = ∂U / ∂P_j.**   (3)

To find a deflection where no load acts, put a **dummy load** Q there, compute U with it, differentiate, then set Q = 0.

**Example: cantilever with a tip load.** M(x) = P(L − x), U = ∫₀ᴸ P²(L − x)²/(2EI) dx = P²L³/(6EI), δ = ∂U/∂P = **PL³/(3EI)**, the result `P05.12` got by integrating twice.

### 2.3 The unit-load method and the flexibility matrix

Differentiate (2) under the integral sign: ∂U/∂P_j = ∫ M (∂M/∂P_j)/EI dx. For a linear structure ∂M/∂P_j = m_j(x), the bending moment caused by a **unit load** at j. So

**δ_j = ∫ M m_j / EI dx  (+ Σ N n_j L/EA + ∫ T t_j / GJ dx),**   (4)

with n_j and t_j the axial force and torque from the unit load. This is the unit-load (virtual-force) method; it is Castigliano with the differentiation done once and for all. Writing the real actions as sums of unit-load actions, M = Σ_k P_k m_k, and so on, gives

**δ = F P,  F_jk = Σ_members ∫ (n_j n_k / EA + m_j m_k / EI + t_j t_k / GJ) dx,**   (5)

and U = ½PᵀFP. F is the **flexibility matrix**; it is plainly symmetric, F_jk = F_kj: the deflection at j due to a unit load at k equals the deflection at k due to a unit load at j (**Maxwell's reciprocal theorem**). It is positive definite because U > 0 for any non-zero load.

### 2.4 Example: the spar under elliptic lift

A cantilever spar of semi-span L carries lift q(x) = q₀√(1 − (x/L)²), the elliptic distribution of an efficient wing, with total lift W = πq₀L/4. The moment at x is M(x) = ∫ₓᴸ q(ξ)(ξ − x) dξ; at the root it is **4WL/(3π)**, the lift acting at its centroid 4L/(3π) from the root. With L = 5 m and W = 18.6 kN: **39.47 kN·m**. For the tip deflection the unit load at the tip gives m = L − x, and (4) is δ = ∫₀ᴸ M(x)(L − x)/EI dx. With EI = 73.1 GPa × 4 × 10⁻⁶ m⁴ = 292.4 kN·m²: **δ_tip = 769 mm**. The same 18.6 kN as a tip load would give WL³/(3EI) = 2651 mm, uniformly spread WL³/(8EI) = 994 mm: the elliptic load is concentrated inboard, so it deflects the tip least. The unbraced number also says that this spar, light enough to be braced, is far too flexible to be a cantilever: 769 mm is 15 % of the span, outside small-deflection theory, and the caps would yield. That is the job of the strut.

### 2.5 Example: an L-bracket in bending and torsion

A bracket built in at its root has a leg of length b, then turns 90° in the same horizontal plane into an arm of length a; a vertical load P acts at the arm's tip. (Think of an offset pitot-boom support, or a gear axle carried on an offset arm.) The arm, with s measured from the tip, bends with M = Ps. The leg, with s measured from the corner, bends with M = Ps and *twists* with T = Pa, because the load acts at a lever arm a from the leg's axis. Energy:

U = ∫₀ᵃ P²s²/(2EI) ds + ∫₀ᵇ P²s²/(2EI) ds + (Pa)²b/(2GJ),

**δ = ∂U/∂P = Pa³/(3EI) + Pb³/(3EI) + Pa²b/(GJ).**   (6)

For a 4130 tube 25.4 × 1.24 mm (I = 6885 mm⁴, EI = 1411 N·m², GJ = 1102 N·m² with G = 80 GPa), a = 0.3 m, b = 0.5 m, P = 500 N: **3.19 + 14.76 + 20.42 = 38.37 mm**. More than half of the deflection is torsion of the leg, a term a "beams only" analysis would miss entirely. With a second load P₂ at the corner, F₁₂ = F₂₁ = b³/(3EI) = 29.52 mm/kN, as Maxwell requires.

### 2.6 Statically indeterminate structures: the redundant from ∂U/∂X = 0

A structure with more supports or members than equilibrium needs is **statically indeterminate**; the extra forces are **redundants**. Release one (cut a member or remove a support) and replace it by an unknown force X. The released structure is determinate, its energy U(X) can be written, and the real X is the one for which the displacement at the cut is zero, by (3):

**∂U / ∂X = 0.**   (7)

That is compatibility, obtained from energy. (It is also the principle of minimum complementary energy: the real force distribution makes U stationary among all those in equilibrium with the loads.)

**Example: the propped cantilever.** A cantilever under uniform load q with a prop at the tip. Take the prop force R as the redundant: M(x) = q(L − x)²/2 − R(L − x). Then ∂U/∂R = ∫ M(−(L − x))/EI dx = 0 gives qL⁴/8 − RL³/3 = 0, **R = 3qL/8**.

### 2.7 The strut-braced spar

Now the engineering problem. The spar is built in at x = 0; the strut runs from a fuselage point a height h below the root to the spar at x = a, at angle θ with sin θ = h/L_s, cos θ = a/L_s, L_s = √(a² + h²). Take the strut tension S as the redundant. The strut pulls the spar down by S sin θ and inboard by S cos θ, so

M(x) = M_q(x) − S sin θ (a − x) for x < a (M_q beyond),  N_spar = −S cos θ on 0 < x < a,  N_strut = S.

The energy has three parts: spar bending, spar axial compression inboard of the strut, and strut tension. (7) gives

**S = sin θ ∫₀ᵃ M_q (a − x)/EI dx / ( sin²θ ∫₀ᵃ (a − x)²/EI dx + cos²θ a/EA_spar + L_s/EA_strut ).**   (8)

Read it as a ratio of flexibilities: the numerator is the spar's deflection at the strut point under the lift (times sin θ), the denominator the total flexibility of the path the strut force takes. A soft strut (EA_strut → 0) carries nothing; a rigid strut on an axially rigid spar gives the propped cantilever, and with a = L and uniform q, S sin θ = 3qL/8.

For the wing of the engineering question: a = 3.0 m, h = 1.2 m (L_s = 3.23 m, sin θ = 0.371), spar EA = 73.1 GPa × 1200 mm², strut EA = 73.1 GPa × 300 mm². Then **S = 31.5 kN** (105 MPa in the strut), whose vertical component **R = 11.7 kN** is 63 % of the lift; the spar's inboard segment is compressed by S cos θ = 29.3 kN (24 MPa). The root moment falls to **4.36 kN·m** and the tip deflection to **48.6 mm**. Compatibility: the spar deflects 15.0 mm at the strut point, and 15.0 mm × sin θ equals the strut's stretch plus the spar's shortening resolved along the strut, to 10⁻⁶. Of the 115 J of strain energy, 73 J is in the strut, 27 J in spar bending and 15 J in spar compression: the strut is not a minor member.

### 2.8 Where the models stop

- Small deflections: the strut angle and the lift direction are assumed fixed. At 49 mm of tip deflection on 5 m that is fine; for the unbraced 769 mm it is not.
- The strut is in tension at positive load factors. At negative g it goes into compression and must be checked for buckling (`P05.15`), usually with a jury strut to shorten its effective length.
- Shear deformation of the spar web, joint flexibility at the strut fittings and the wing's torsion are left out; each adds flexibility and changes S by a few per cent.
- Castigliano's second theorem needs a linear structure. With plasticity or large deflections U is not quadratic in the loads, and (3) must be replaced by the more general complementary-energy form.

## 3. Build the computational model

Open `starter/main.py`. The `Member` dataclass, `elliptic_lift`, `truss_lengths` and `main` are given.

### Step 1: energy and Castigliano by differentiation

`strain_energy_axial`, `strain_energy_torsion` are (1) and its torsion twin; `strain_energy_bending(x, M, EI)` is (2) with `scipy.integrate.simpson`. `castigliano_deflection(energy, loads, j)` is a central difference of a user-supplied U(P) with respect to P_j. Because U is quadratic, the central difference has no truncation error, only round-off.

### Step 2: unit loads and the flexibility matrix

`unit_load_deflection(x, M, m, EI)` is (4) for bending. A `Member` holds L, EA, EI, GJ and a function `actions(s)` returning, for each unit load, the axial force, moment and torque along the member, shape (n_loads, 3, n_stations). `flexibility_matrix(members)` is (5):

```python
a = mb.actions(s)                                       # (n_loads, 3, n)
inv = np.array([1.0 / mb.EA, 1.0 / mb.EI, 1.0 / mb.GJ])  # 1/inf = 0 drops a term
Fm = simpson(np.einsum("jcs,kcs,c->jks", a, a, inv), x=s, axis=-1)
```

`strain_energy_members` integrates the real N²/2EA + M²/2EI + T²/2GJ directly, which the tests compare with ½PᵀFP. `l_bracket_members` builds the two members of 2.5; `l_bracket_tip_deflection` is (6).

### Step 3: the spar

`cantilever_moment(x, q)` computes M(x) = S₁(x) − xS₀(x) from cumulative integrals of q and xq from the tip (cumulative trapezoid). `deflection_by_integration(x, M, EI)` integrates M/EI twice from the root with w(0) = w'(0) = 0: the direct route of `P05.12`, and the independent check of the unit-load result.

### Step 4: trusses, determinate and not

`truss_equilibrium_matrix` builds C with C N = F at the free degrees of freedom, N positive in tension (the method of joints of `P05.03` as one matrix: at each joint the bars' pulls toward their far ends balance the load, so a column holds minus the bar's unit vectors). Its transpose maps joint displacements to bar elongations, which is what makes the stiffness matrix C diag(EA/L) Cᵀ. `truss_force_method` finds a particular solution N₀ by least squares, the self-stress states S as the null space of C (from the SVD), and the redundants X = −(SᵀDS)⁻¹SᵀDN₀ with D = diag(L/EA), which is (7) for every redundant at once. `truss_unit_load_displacement` is (4) for bars, with the unit-load forces from any equilibrium solution (the virtual-force principle needs only static admissibility, so even for an indeterminate truss a particular solution will do). `truss_stiffness_method` is the independent check: K = C diag(EA/L) Cᵀ, u = K⁻¹F.

### Step 5: the braced spar

`braced_spar(L, a, h, EI, EA_spar, EA_strut, q)` evaluates (8) by Simpson's rule, builds the final M(x), integrates the deflection, and reports the strut force, the reactions, the tip deflection and the compatibility residual.

Running `main()` prints

```
cantilever spar L = 5.0 m, EI = 292.4 kN m^2, elliptic lift 18.6 kN: root moment 39.47 kN m (4 W L / 3 pi = 39.47); tip deflection unit load 768.95 mm, direct integration 768.95 mm
  same total as a tip load: 2650.5 mm; uniform: 993.9 mm

L-bracket a = 0.3 m, b = 0.5 m, P = 500 N: closed form 38.373 mm, flexibility 38.373 mm, dU/dP 38.373 mm
  arm bending 3.188 mm, leg bending 14.760 mm, leg torsion 20.424 mm; F12 = 29.5202, F21 = 29.5202 mm/kN

engine mount, 2 bars (determinate): bar forces force method -4.193, +7.019 kN; stiffness method -4.193, +7.019 kN
  vertical deflection: unit load -0.6291 mm, stiffness -0.6291 mm

engine mount, 3 bars (1 redundant): bar forces force method -4.994, +1.441, +6.144 kN; stiffness method -4.994, +1.441, +6.144 kN
  vertical deflection: unit load -0.6170 mm, stiffness -0.6170 mm

braced spar, strut at a = 3.0 m, h = 1.2 m: strut tension 31.51 kN, vertical part 11.70 kN; root moment 4.36 kN m (unbraced 39.47); tip 48.6 mm (unbraced 768.9); compatibility residual 9.1e-07
  rigid strut at the tip, uniform 1 kN/m: R = 1.8750 kN against 3 q L / 8 = 1.8750 kN
```

followed by the two sweep tables of part 4, and writes `energy_methods.png` (bending moment and deflection of the unbraced and braced spar, with the unit-load tip deflection marked).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-16-energy-methods`. The tests:

1. **Energy formulas** (1), (2) and the torsion form; ∫M²/2EI for a tip load equals P²L³/6EI.
2. **Castigliano**: the numerical ∂U/∂P equals PL³/3EI to 10⁻⁹ at two load levels.
3. **Unit load against direct integration** for elliptic, uniform and linearly varying lift, to 10⁻⁵; M(0) equals the moment of the lift about the root and M(L) = 0 (equilibrium).
4. **Closed forms**: qL⁴/8EI and q(L − x)²/2 for uniform load; total elliptic lift and 4WL/3π.
5. **L-bracket**: the flexibility matrix and Castigliano equal (6); the torsionally rigid limit leaves (a³ + b³)/3EI.
6. **Flexibility matrix**: F₁₂ = F₂₁ = b³/3EI, positive definite, and U = ½PᵀFP on random loads.
7. **Superposition** of deflections.
8. **Trusses** (determinate and one redundant): force method = stiffness method to 10⁻⁹; C N = F (joint equilibrium) to 10⁻⁸; unit-load displacements = stiffness displacements.
9. **Redundant limit**: a bar of vanishing EA carries nothing and the other two revert to the determinate solution; the extra bar stiffens the truss.
10. **Braced spar**: compatibility residual below 10⁻⁴; root moment and shear from global equilibrium; M(L) = 0; strut in tension.
11. **Rigid limit**: zero deflection at a rigid strut; R = 3qL/8 for a = L.
12. **Monotonicity**: S rises from 0 with strut stiffness; the tip deflection falls.
13. **main** writes the figure.

### Compare

The direct stiffness method is the established tool for trusses and frames (it is what a finite-element code does), and it agrees with the force method to round-off on both trusses. The difference between the two is the choice of unknowns: the stiffness method solves for displacements, one per free degree of freedom (two here); the force method for redundant forces, one per degree of indeterminacy (zero or one here). For a structure with many joints and few redundants, the force method is the smaller problem, which is why it was the method of choice for aircraft structures before computers, and why it remains the method for hand checks of an FE model.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Unit load against direct integration.** Tip deflection of the 5 m spar under elliptic, uniform and tip-concentrated lift of the same total, by (4) and by double integration. Predict the ordering of the three deflections and the agreement between methods.
2. **The L-bracket, term by term.** The three terms of (6) for a = 0.3 m and b from 0.1 to 0.8 m. Predict at which b the torsion term overtakes the leg's bending term (set Pa²b/GJ = Pb³/3EI).
3. **The truss with and without a redundant bar.** The engine-mount truss of `main` (nodes at (0, 0), (0, 0.3), (0, 0.6) and (0.5, 0.25) m, load 2.0 kN forward and 5.9 kN down at the engine node, 4130 tubes 19.05 × 0.889 mm) with two bars and with a third, middle bar. Predict whether the third bar is in tension or compression and how much it reduces the deflection. (Bar forces are printed in the order bottom, middle, top; tension positive.)
4. **Strut stiffness sweep.** EA_strut from 0 to 10⁶ times the design value. Predict the limits of S, the root moment and the tip deflection.
5. **Strut position sweep.** a from 1.5 to 4.0 m at fixed h. Predict where the largest bending moment anywhere on the spar is smallest.

Raw data format: CSV tables `load_case, delta_unit_load_mm, delta_integration_mm` (1); `b_m, arm_mm, leg_bending_mm, leg_torsion_mm` (2); `case, N1_kN, N2_kN, N3_kN, uy_mm` (3); `EA_ratio, S_kN, R_kN, Mroot_kNm, tip_mm` (4); `a_m, S_kN, R_kN, Mroot_kNm, Mmax_kNm, tip_mm` (5); plus `energy_methods.png`.

What you should find. (1) 769, 994 and 2651 mm, the two methods agreeing to six figures. (2) The torsion term equals the leg bending when b² = 3a²EI/GJ, b = 0.59 m here; beyond it torsion dominates. (3) The middle bar takes 1.44 kN of tension, the top bar's tension falls from 7.02 to 6.14 kN and the bottom bar's compression rises from 4.19 to 4.99 kN; the vertical deflection falls by 1.9 %, from 0.629 to 0.617 mm. A redundant bar mostly redistributes load; here it is nearly horizontal while the load is mostly vertical, so it adds little vertical stiffness. (4) and (5):

```
EA_strut / EA_ref   S [kN]    R [kN]   root M [kN m]   tip [mm]
                 0     0.00     0.00           39.47     768.9
              0.01     7.33     2.72           31.30     601.3
               0.1    24.24     9.00           12.46     214.7
                 1    31.51    11.70            4.36      48.6
                10    32.48    12.06            3.28      26.3
             1e+06    32.59    12.11            3.15      23.7

strut position a [m]   S [kN]   R [kN]   root M [kN m]   max |M| [kN m]   tip [mm]
                 1.5    46.67    29.15           -4.26           16.86     264.2
                 2.0    39.72    20.44           -1.41           11.61     163.2
                 2.5    35.05    15.17            1.55            7.45      93.4
                 3.0    31.51    11.70            4.36            4.36      48.6
                 3.5    28.63     9.28            6.98            6.98      23.5
                 4.0    26.19     7.53            9.36            9.36      13.9
```

The design strut already has 97 % of the rigid-strut force: beyond a modest stiffness, making the strut stiffer buys little. Moving the strut outboard shortens the overhang and the tip deflection keeps falling, but the root moment rises; inboard, the overhang's moment grows. The largest moment anywhere on the spar is smallest near a = 3.0 m, where the root and strut-point moments are about equal (4.36 and 4.32 kN·m). That is the classic placement of a lift strut a little beyond mid-span. (At a = 1.5 and 2.0 m the root moment is negative: the strut lifts the inboard spar so hard that it bends the other way.)

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| root moment, elliptic lift | 4WL/3π = 39.47 kN·m | 39.47 kN·m | < 10⁻³ | √ singularity of q at the tip, trapezoid |
| tip deflection, unit load vs integration | equal | 768.95 / 768.95 mm | < 10⁻⁵ | quadrature |
| Castigliano vs PL³/3EI | equal | 10⁻⁹ | — | round-off |
| L-bracket | (6) = 38.373 mm | flexibility and ∂U/∂P 38.373 mm | < 10⁻¹⁰ | — |
| Maxwell | F₁₂ = F₂₁ = b³/3EI | 29.5202 / 29.5202 mm/kN | < 10⁻¹² | — |
| truss bar forces | stiffness method | force method | < 10⁻⁹ | — |
| propped cantilever | 3qL/8 = 1.8750 kN | 1.8750 kN | < 10⁻⁶ | — |
| braced spar compatibility | w(a) sin θ = stretch + shortening | residual 9 × 10⁻⁷ | — | Simpson against trapezoid |

Every check passes at the level of its quadrature. The checks are independent in the sense that matters: direct integration never uses the unit load, the stiffness method never uses the force method's null space, and compatibility is checked with a deflection computed after S was fixed by energy. The engineering answer is **strut tension 31.5 kN (vertical 11.7 kN, 63 % of the lift), root moment 4.4 kN·m against 39.5 unbraced, tip deflection 49 mm**, on the assumptions of 2.8.

What to improve: the strut-fitting flexibility (a few per cent on S), the negative-g case with the strut in compression (a buckling check), and a spanwise-varying EI for a tapered spar (`braced_spar` already accepts an array EI if you pass one). What not to improve: the number of stations, already ten times finer than the uncertainty in the lift distribution.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (6) and of 3qL/8; `energy_methods.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a strut position for the wing, with the strut force, the spar's largest moment and the tip deflection, and the assumption most likely to change the answer.

### Adversarial self-checks

- **What if the lift doubles?** Everything doubles: S, M, the deflections. The energy quadruples. Nothing changes shape, because the structure is linear; the first thing to change at higher load is that the spar caps yield or the strut goes into compression at negative g.
- **What if the strut stiffness goes to zero, or to infinity?** Zero: S = 0 and the unbraced cantilever returns. Infinity: S tends to 32.6 kN, not to the propped-cantilever value, because the spar is still axially flexible inboard of the strut and the strut is not at the tip. Both limits are in the sweep.
- **Which assumption fails first in the real world?** The pinned strut end and the built-in spar root: real fittings have bolt clearances and flexibility, which act as springs in series with the strut and lower S. Then the fixed strut angle, if deflections grow.
- **What measurement could falsify the model?** A strain gauge on the strut and a dial gauge at the tip during a ground proof-load test with sandbags (how such wings are tested): S and the tip deflection against (8) and the computed deflection.
- **Does the result survive a change of method?** The force method and the stiffness method agree to 10⁻⁹ on the truss; the unit load and direct integration agree to 10⁻⁵ on the spar; the strut force from (8) satisfies compatibility to 10⁻⁶.
- **Which decision changes if the spar EI is 20 % lower than assumed?** The spar becomes more flexible relative to the strut, so the strut takes more load (S rises toward its rigid limit) and the root moment falls slightly; the tip deflection rises. The strut, not the spar, would need the extra margin.

## Check yourself

1. Write the strain energy of a bar in tension, a shaft in torsion and a beam in bending.
2. State Castigliano's second theorem. How do you find a deflection at a point where no load acts?
3. What is the unit-load integral for the deflection of a beam, and what is m(x)?
4. State Maxwell's reciprocal theorem in words, and say which property of the flexibility matrix expresses it.
5. A cantilever under uniform load q has a rigid prop at its tip. What is the prop force, and which equation gives it in the energy method?
6. Why does an L-bracket loaded at its tip deflect more than the sum of the bending deflections of its two legs?
7. In the braced spar, what happens to the strut force as the strut becomes very soft, and why?
8. What does the force method solve for, and when is it a smaller problem than the stiffness method?
