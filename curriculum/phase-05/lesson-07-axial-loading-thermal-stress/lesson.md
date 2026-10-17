# P05.07  Axial Loading & Thermal Stress

> A bolt is a spring stretched against a spring squeezed. Load them, cool them, or tighten them carelessly, and the joint tells you which of the three you forgot.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the elongation of uniform, stepped and tapered bars, and solve statically indeterminate bars by adding a compatibility condition to equilibrium, with mechanical and thermal loads together. 2. Implement a 1-D bar finite-element solver with initial strains and use it as the independent check of every closed form in the lesson, including second-order convergence on a tapered bar. 3. Model a preloaded bolted joint: bolt stiffness (shank and thread in series), member stiffness by the pressure-cone model integrated in closed form, the joint diagram, load sharing, separation and the change of preload with temperature. 4. Decide whether a joint meets its requirements over a temperature range when the preload itself is uncertain, using a seeded Monte Carlo study of torque–tension scatter. |
| **Prerequisites** | `P05.06` (stress, strain, Hooke's law, the material table), `P05.02` (equilibrium, determinacy), `P05.04` (screw-thread friction, the origin of the torque–tension relation), `P03.13` (the 1-D bar finite element), `P03.04` (quadrature), `P01.13` (linear systems), `P03.17` and `P02.16` (Monte Carlo and random variables). Skill check: you can write the stiffness of two springs in series and in parallel. |
| **Engineering question** | An M8 × 1.25 Ti-6Al-4V bolt clamps two 10 mm 7075-T6 fittings of an engine-pylon attachment (an illustrative joint) with a preload of 20 kN. Each bolt sees a limit separating load of 12 kN (18 kN ultimate). The joint is assembled at 20 °C, soaks to −55 °C in cruise and to +70 °C on a hot ramp. Does it stay clamped at ultimate load when cold, does the bolt stay below yield at limit load when hot, and does the answer survive the scatter of torque tightening? The answers: the bolt takes C = 0.308 of any external load; cold, the preload falls to 16.0 kN and the joint separates at 23.1 kN, above 18 kN; hot, the bolt stress at limit is 720 MPa, below 880 MPa; but with a nut-factor scatter of 12.5 % (one standard deviation), 4.2 % of joints separate below ultimate when cold and 3.6 % exceed yield at limit when hot. Torque control is not good enough for this joint. |
| **Assumptions and limits** | Linear elasticity (`P05.06`); uniaxial stress in bars (no stress concentration at shoulders and threads); uniform temperature in each part, with properties independent of temperature; the bolt and members as two linear springs; the members' stiffness from a 30° pressure cone, an engineering model that is reasonable when the members extend well beyond the cone and poor near an edge or with several bolts close together; the external load introduced at the head and nut faces (the simplest and most conservative loading plane for the bolt); no bending of the bolt, no prying, no gasket, no relaxation or embedding. Material values are typical, not design allowables. |
| **Analytical work** | δ = NL/EA from Hooke's law and the method of sections; stepped bars as springs in series; the tapered bar by integration; the constrained bar −EαΔT; the fixed–fixed two-segment bar by compatibility; the tensile stress area and bolt stiffness; the frustum stiffness by partial fractions; the joint diagram from equilibrium and compatibility; the separation load; the thermal change of preload; the deterministic preload window that meets both cold and hot requirements. |
| **Computational deliverable** | `starter/main.py`: `axial_deformation`, `stepped_bar_deformation`, `tapered_bar_deformation`, `constrained_thermal_stress`, `bar_fe_solve`, `fixed_fixed_two_segment`, `tensile_stress_area`, `bolt_stiffness`, `_cone_compliance`, `frustum_stiffness`, `member_stiffness`, `joint_constant`, `separation_load`, `joint_loads`, `thermal_preload_change`, `preload_from_torque` (`MATERIALS`, the `Joint` dataclass, the lesson joint and `main` are given). Tests prove: the FE solver satisfies equilibrium and reproduces the stepped bar, the fixed–fixed bar with load and temperature, and the constrained thermal stress to round-off, and the tapered bar at order 2.0 ± 0.1; superposition of load and temperature; the frustum stiffness against `scipy.integrate.quad` to 10⁻¹⁰ and its thin-plate limit; layered member stiffness against an independent quadrature of the double cone to 10⁻⁸; the joint loads satisfy F_b − F_m = P and (F_b − F_i)/k_b = (F_i − F_m)/k_m and match a two-element FE model of the joint with a misfit; separation and the rigid-member and rigid-bolt limits; the thermal preload change vanishes for matched materials, is linear in ΔT and matches FE. |
| **Practical activity** | **computational experiment**: FE convergence on a tapered bar, a joint parameter study (bolt material, temperature, load), and a seeded Monte Carlo study of torque–tension scatter against the cold-separation and hot-yield requirements. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `bolted_joint.png`, the derivations of the fixed–fixed bar, the frustum stiffness, the joint diagram and the thermal preload change, the experiment tables, and a half-page note answering the engineering question with a recommended preload method. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the joint diagram is derived from equilibrium and compatibility (not quoted); all public tests pass; the FE check is reported for every closed form; the note states the cold separation load, the hot bolt stress and the Monte Carlo fractions with their assumptions. |
| **References** | `@openstax-university-physics-vol2` 1.3 Thermal Expansion: "Linear Thermal Expansion" and "Thermal Stress" (required-free): α and the restrained-bar stress. `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Trusses": "Deflections" and "Statically indeterminate trusses" (required-free): δ = PL/AE and compatibility for redundant bars. `@nasa-rp-1228-fastener-design-manual` "General Design Information" → "Fastener Torque"; "Fatigue-Resistant Bolts" → "Fatigue (Cyclic) Loading of Bolts"; "Design Criteria" (required-free, https://ntrs.nasa.gov/citations/19900009424): the torque coefficient relation, bolt and member spring constants and the load–elongation diagram with preload, external load and separation. `@nasa-tm-x-73305-astronautic-structures-manual` Vol. III, Section D "THERMAL STRESSES": 3.0.1 Unrestrained Beam -- Thermal Loads Only; 3.0.2 Restrained Beam -- Thermal Loads Only (supplementary-free, https://ntrs.nasa.gov/citations/19760071127). `@tudelft-alderliesten-aerospace-structures-materials` 2.2 The effect of ambient temperature; 11.2 Mechanically fastened joints (required-free, link only: https://books.open.tudelft.nl/home/catalog/book/197). |
| **Extensions** | ★ Make the moduli and expansion coefficients functions of temperature (aluminium loses several per cent of its modulus by 150 °C) and redo the hot case. ★ Move the loading plane inside the members (a load introduction factor n < 1) and show that the bolt's share falls to nC. ★ Replace the pressure cone by a 2-D axisymmetric FE model of the members from `P03.13` and find the cone angle that matches it. |
| **Time** | unmeasured |

## 1. The engineering problem

An aircraft is held together by tens of thousands of fasteners, and the critical ones (wing-to-fuselage, engine pylon, landing-gear fittings) are preloaded bolts. A preloaded bolt does not work the way intuition suggests. Tighten it to 20 kN, then hang 12 kN on the joint, and the bolt tension does not become 32 kN or stay at 20 kN: it rises to 23.7 kN, because most of the external load goes into *unloading* the clamped plates rather than stretching the bolt further. That is why preload protects a bolt against fatigue: the alternating part of the load the bolt feels is only a fraction of the alternating load on the joint. It works until the plates lose contact. After separation the bolt carries everything, the joint rattles, and fatigue cracks follow.

Temperature changes the preload without anyone touching the nut. Titanium bolts in aluminium fittings are common because titanium is strong and light; but aluminium expands almost three times as much per kelvin. Cool the joint by 75 K in cruise and the plates shrink more than the bolt, the bolt relaxes, and a fifth of the preload disappears. Heat it on a ramp in the sun and the preload rises. On top of that, the preload is set by turning a wrench to a torque, and the torque–tension relation depends on friction in the threads and under the head, which scatters by tens of per cent.

This lesson builds the tools to answer that question with numbers. Part 2 starts with the axial bar (δ = NL/EA, the constrained bar, the indeterminate bar) because a bolted joint *is* two bars, then derives the joint diagram. Part 3 implements a 1-D bar finite-element solver from `P03.13` that checks every closed form independently. Part 4 sweeps the joint over temperature, bolt material and preload scatter. Success criteria, written first: every closed form equals the FE result to round-off (or converges to it at second order); the joint loads satisfy equilibrium and compatibility to round-off; and the note states the probability that the joint violates each requirement under a stated model of torque scatter.

The organising idea of the phase is *equilibrium, compatibility and constitutive law*, and this lesson needs all three in every result: equilibrium of a cut bar or a joint, compatibility of the deformations of parts that are connected, and Hooke's law with thermal strain, ε = σ/E + αΔT, from `P05.06`.

## 2. Understand and derive

### 2.1 The axially loaded bar

Cut a bar of length L, area A and modulus E anywhere along its length; equilibrium of either piece gives the internal force N (tension positive). The stress is N/A, the strain (Hooke's law) N/(EA), and the elongation, adding up the strain along the bar (compatibility),

**δ = N L / (E A).**   (1)

The bar is a spring of stiffness **k = EA/L**. A 2024-T3 rod 1.2 m long, 100 mm² in area, carrying 15 kN stretches 2.462 mm. Segments in series carrying their own internal forces add their elongations, δ = Σ NᵢLᵢ/(EᵢAᵢ), which is the stepped bar; when the area varies continuously, the sum becomes an integral. For an area that falls linearly from A₀ to A₁ over L, carrying P,

**δ = ∫₀ᴸ P dx / (E A(x)) = P L ln(A₀/A₁) / (E (A₀ − A₁)),**   (2)

which tends to PL/(EA₀) as A₁ → A₀. A 2024-T3 bar tapering from 200 to 50 mm² over 500 mm under 20 kN stretches 1.26429 mm. Equation (1) is Saint-Venant's idealisation: near a load introduction, a shoulder or a thread the stress is not uniform, and (1) is only correct a couple of widths away.

### 2.2 Thermal strain and the restrained bar

A free bar heated by ΔT grows by αLΔT, with α the coefficient of linear expansion (OpenStax 1.3): about 23 × 10⁻⁶ K⁻¹ for aluminium alloys, 12 × 10⁻⁶ K⁻¹ for steel and 8.6 × 10⁻⁶ K⁻¹ for Ti-6Al-4V. Hooke's law with a thermal strain is

**ε = σ/E + α ΔT.**   (3)

If the bar is held between rigid walls, compatibility says ε = 0, and (3) gives

**σ = −E α ΔT,**   (4)

compressive when heated, independent of the length and the area. Heating by 50 K gives −84.8 MPa in 2024-T3, −125.0 MPa in 4130 steel and −48.9 MPa in Ti-6Al-4V: steel's higher modulus outweighs its smaller expansion. Equation (4) is the reason cryogenic propellant lines have bellows (Problem 3) and the reason a bar between two parts at different temperatures is never stress-free.

### 2.3 Statically indeterminate bars

Join a 4130 steel segment (300 mm, 80 mm²) to a 2024-T3 segment (400 mm, 150 mm²) and fix both outer ends to rigid walls. Push the junction with 10 kN along the bar. Equilibrium of the junction gives one equation, N₂ − N₁ + P = 0, for two unknown internal forces: the bar is statically indeterminate (`P05.02`). The missing equation is **compatibility**: the two segments together cannot change length. With the junction displacement u and the segment stiffnesses kᵢ = EᵢAᵢ/Lᵢ, segment 1 stretches by u and segment 2 by −u, each by a mechanical plus a thermal part (3):

u = N₁/k₁ + α₁L₁ΔT,   −u = N₂/k₂ + α₂L₂ΔT.

Substituting into equilibrium,

**u = (P + k₁α₁L₁ΔT − k₂α₂L₂ΔT) / (k₁ + k₂),  N₁ = k₁(u − α₁L₁ΔT),  N₂ = k₂(−u − α₂L₂ΔT).**   (5)

With ΔT = 0: u = 0.1218 mm, N₁ = 6660 N (83.3 MPa tension) and N₂ = −3340 N (22.3 MPa compression). The load splits in proportion to stiffness, not area or length alone. Heating both segments by 40 K as well: u = 0.0954 mm, N₁ = −2790 N (−34.9 MPa), N₂ = −12 790 N (−85.3 MPa). The solution is linear, so load and temperature superpose (the tests check it).

### 2.4 The 1-D bar finite element, with initial strain

`P03.13` assembled a stiffness matrix from bar elements, each (EA/L)[[1, −1], [−1, 1]]. An initial strain ε₀ (thermal αΔT, or a misfit: a part made too short) enters as an equivalent pair of nodal forces EAε₀(−1, +1), and the element force is N = EA((u_j − u_i)/L − ε₀). Solving K u = f + f₀ with the fixed displacements removed gives the displacements; the reactions are R = K u − f₀ − f at the fixed nodes. For bars of constant section with nodal loads the method is exact, so it must reproduce (5) to round-off; for a tapered bar with midpoint areas it converges at second order to (2). This solver is the independent check of every closed form in this lesson, including the joint.

### 2.5 The bolt and the members as springs

**The bolt.** Within the grip (the clamped length l), the bolt has an unthreaded shank of area A_d = πd²/4 and length l_s, and usually some thread of length l_t, whose effective area is the **tensile stress area**, for ISO metric threads

**A_t = (π/4)(d − 0.9382 p)²**   (6)

(p the pitch): 36.61 mm² for M8 × 1.25, against 50.27 mm² for the shank. The two parts are springs in series, 1/k_b = l_s/(EA_d) + l_t/(EA_t). The lesson bolt (Ti-6Al-4V, l_s = 14 mm, l_t = 6 mm) has k_b = 257.2 kN/mm.

**The members.** The clamping force spreads out from under the head and the nut into the plates. The pressure-cone model replaces the compressed region by two cones (frusta) of half-angle α = 30°, starting at the bearing-face diameter D_w and meeting at the mid-plane of the grip. At distance s from a face the cone's diameter is D(s) = D_w + 2s tan α and the loaded area is an annulus π(D² − d²)/4 around the hole. A frustum of height t is a tapered bar, and its compliance is

∫₀ᵗ ds / (E π((D(s)/2)² − (d/2)²)).

Partial fractions, 1/(u² − a²) = (1/2a)[1/(u − a) − 1/(u + a)] with u = D(s)/2 and a = d/2, give the closed form

**k_frustum = π E d tan α / ln[ (2t tan α + D_w − d)(D_w + d) / ((2t tan α + D_w + d)(D_w − d)) ].**   (7)

For the lesson joint (7075-T6, d = 8 mm, D_w = 12 mm, two 10 mm plates), each frustum has k = 1153.7 kN/mm and the two in series give **k_m = 576.8 kN/mm**, more than twice the bolt. For a stack of different materials, each layer on each side of the mid-plane is a slice of a frustum, with the same antiderivative evaluated between its limits, and all slices add as compliances. The cone angle is a model choice, not a law: the angle that best represents a joint depends on its geometry, and changing it by ±5° moves k_m to 511.9 or 648.0 kN/mm (about ±12 %) and C to 0.334 or 0.284. The model fails when the cone would extend past the edge of the part or overlap the next bolt's cone.

### 2.6 The joint diagram

Tighten the nut. The bolt stretches by F_i/k_b and the members shorten by F_i/k_m; both carry the preload F_i, the bolt in tension and the members in compression (equilibrium of the head). Now apply an external separating load P at the head and nut faces. Equilibrium of the bolt head with the plate surface:

**F_b − F_m = P,**   (8)

with F_b the bolt tension and F_m the clamping force (compression, as a positive number). Compatibility: the head moves outward by the same amount whether you follow the bolt or the members, so the bolt's extra stretch equals the members' recovery,

**(F_b − F_i)/k_b = (F_i − F_m)/k_m.**   (9)

Solving (8) and (9):

**F_b = F_i + C P,   F_m = F_i − (1 − C) P,   C = k_b/(k_b + k_m).**   (10)

C is the **joint constant**, the bolt's share of the external load: 0.308 for the lesson joint. The bolt feels less than a third of every load change; the members' unloading takes the rest. The clamping force reaches zero at the **separation load**

**P_sep = F_i / (1 − C),**   (11)

28.92 kN at assembly. Beyond it the bolt carries all of P alone. Plotting the bolt's line (slope k_b) and the members' line (slope −k_m, starting from the total misfit) gives the classic joint diagram; `bolted_joint.png` draws it. Two limits check (10): rigid members (k_m → ∞) give C = 0, the bolt feels nothing until separation; a rigid bolt gives C = 1, it takes everything.

For the lesson joint at 20 °C and the 12 kN limit load: F_b = 23.70 kN (647 MPa on the tensile stress area) and F_m = 11.70 kN. At the 18 kN ultimate load: F_b = 25.55 kN (698 MPa), F_m = 7.55 kN.

### 2.7 Temperature and preload

When the joint changes temperature by ΔT, the members try to grow by Σαᵢtᵢ ΔT and the bolt by α_b l ΔT. The difference, the thermal misfit δ_th = (Σαᵢtᵢ − α_b l)ΔT, must be taken up elastically by the bolt and members in series, so

**ΔF_i = (Σ αᵢ tᵢ − α_b l) ΔT · k_b k_m / (k_b + k_m).**   (12)

For the lesson joint cooled by 75 K: δ_th = (23.6 − 8.6) × 10⁻⁶ × 20 mm × (−75 K) = −22.5 μm, the series stiffness is 177.9 kN/mm, and **ΔF_i = −4.00 kN**: the preload falls from 20.0 to 16.0 kN. Heated by 50 K it rises by 2.67 kN to 22.67 kN. Then, from (10) and (11):

| case | preload (kN) | P_sep (kN) | F_b at limit (kN) / stress (MPa) | F_b at ultimate (kN) / stress (MPa) | clamp left at ultimate (kN) |
|---|---|---|---|---|---|
| assembly 20 °C | 20.00 | 28.92 | 23.70 / 647 | 25.55 / 698 | 7.55 |
| cold −55 °C | 16.00 | 23.13 | 19.70 / 538 | 21.55 / 589 | 3.55 |
| hot +70 °C | 22.67 | 32.78 | 26.37 / 720 | 28.22 / 771 | 10.22 |

With the nominal preload, every requirement is met: the joint stays clamped at ultimate even when cold (23.13 > 18 kN), and the hot bolt stays below the typical 880 MPa yield at limit and below the 950 MPa ultimate at ultimate. A 4130 steel bolt in the same joint (stiffer: k_b = 463.4 kN/mm, C = 0.445) loses even more preload when cold (to 15.61 kN): the smaller mismatch is outweighed by the higher stiffness in (12). (A real steel aerospace bolt is heat-treated far above the typical 4130 strength in the table; the comparison is about stiffness and expansion only.)

### 2.8 Torque and preload

The preload is set by a wrench. The short-form relation

**T = K F_i d**   (13)

lumps thread friction, head friction and thread geometry (`P05.04`) into the nut factor K, about 0.2 for dry steel threads; the Fastener Design Manual's "Fastener Torque" section gives the relation and the ways torque is measured, and its coefficients are typical values, not data for a particular fastener and lubricant. 20 kN on an M8 bolt needs T = 32.0 N·m at K = 0.2. Because friction scatters, the same torque gives different preloads: with K normally distributed about 0.20 with a standard deviation of 0.025 (an assumption of this lesson, ±25 % at two standard deviations), 98 % of joints torqued to 32 N·m have preloads between 15.5 and 28.3 kN.

### 2.9 Where the models stop

- The joint diagram is linear: real joints embed (surface roughness flattens, losing a few per cent of preload in the first hours), and a gasket or sealant adds a soft spring.
- Loading at the head and nut faces is the worst case for the bolt's share; load introduced inside the members lowers it.
- Bending: an eccentric load or a joint that opens on one side (prying) adds bending stress to the bolt, which (10) ignores.
- Properties change with temperature; at +70 °C the change is small for these alloys, at 150 °C it is not.
- The Monte Carlo result is only as good as the assumed scatter of K; the method transfers, the percentages do not.
- The values in `MATERIALS` (E, α, Fty, Ftu for the four alloys) are typical room-temperature values rounded from producers' data sheets, with the status explained in `P05.06`; the expansion coefficients are representative of 20–100 °C.

## 3. Build the computational model

Open `starter/main.py`. `MATERIALS`, the `Joint` dataclass (with `k_b`, `k_m`, `A_t` and `preload_at(dT)`), the lesson joint `LESSON_JOINT`, the load case and `main` are given.

### Step 1: bars

`axial_deformation`, `stepped_bar_deformation`, `tapered_bar_deformation` and `constrained_thermal_stress` are (1), its sum, (2) and (4).

### Step 2: the FE solver

```python
def bar_fe_solve(x, elements, EA, loads, fixed, eps0=None):
    # elements: (m, 2) node pairs, or None for consecutive nodes
    for (i, j), ea, L, e in zip(el, EA, Ls, e0):
        K[np.ix_([i, j], [i, j])] += ea / L * np.array([[1, -1], [-1, 1]])
        f0[i] -= ea * e; f0[j] += ea * e          # initial strain as nodal forces
    u[free] = np.linalg.solve(K[np.ix_(free, free)], (loads + f0)[free])
    N = EA * ((u[el[:, 1]] - u[el[:, 0]]) / Ls - e0)
    R = K @ u - f0 - loads; R[free] = 0.0
```

Allowing any connectivity lets two elements share the same two nodes: that is a bolt and its members in parallel. `fixed_fixed_two_segment` is (5).

### Step 3: the bolt and the members

`tensile_stress_area` is (6); `bolt_stiffness` the two springs in series. `_cone_compliance(s_a, s_b, E, d, D_w, tan_a)` is the antiderivative of 2.5 evaluated between two distances from a face; `frustum_stiffness` calls it from 0 to t; `member_stiffness` walks through the layers from the head side, splitting each at the mid-plane and measuring s from the nearer face.

### Step 4: the joint

`joint_constant`, `separation_load`, `joint_loads` ((10), with F_b = P and F_m = 0 after separation, vectorised over P), `thermal_preload_change` (12) and `preload_from_torque` (13).

Running `main()` prints the bar results of 2.1–2.3 with the FE differences (4.5 × 10⁻¹³ mm and below), the tapered-bar convergence table

```
elements   FE delta (mm)   rel. error
       2   1.20263        4.88e-02
       4   1.24593        1.45e-02
       8   1.25940        3.87e-03
      16   1.26305        9.84e-04
      32   1.26398        2.47e-04
```

the joint table of 2.7, the FE check of the joint, the steel-bolt comparison and the Monte Carlo lines, and writes `bolted_joint.png` (the joint diagram; bolt and clamp forces against external load at the three temperatures).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-07-axial-loading-thermal-stress`. The tests:

1. **Bars**: (1); equal halves in series equal one bar; a compressive segment shortens.
2. **FE against sections**: a three-segment bar with intermediate loads gives the internal forces of the method of sections to 10⁻¹², reactions balancing the loads, and the stepped-bar elongation.
3. **Fixed–fixed bar**: (5) against FE to 10⁻¹² mm for four combinations of load and temperature, with junction equilibrium, overall equilibrium and zero total elongation.
4. **Superposition** of load and temperature in (5).
5. **Constrained thermal stress**: FE on a bar between walls gives −EαΔT in every element, for any length and area; a free bar gives zero stress and αΔTL.
6. **Tapered bar**: (2) against `quad`; FE errors converging at order 1.9–2.1; the uniform limit.
7. **Bolt**: A_t of M8 between 36 and 37 mm²; all-shank stiffness EA/l; threads softer; twice the grip, half the stiffness.
8. **Frustum**: (7) against `quad` of the cone integral to 10⁻¹⁰; a very thin frustum is an annulus.
9. **Members**: symmetric joint = half a frustum's stiffness; splitting a layer changes nothing; a three-material stack against an independent quadrature to 10⁻⁸; proportional to E.
10. **Joint diagram**: (8) and (9) to round-off over a range of P; slope C; separation, continuity at separation, the rigid-member and rigid-bolt limits.
11. **Joint against FE**: two parallel elements with the bolt made short by the misfit F_i(1/k_b + 1/k_m) reproduce F_b and F_m to 10⁻¹⁰ at four loads.
12. **Thermal preload**: zero for matched materials, negative cold and positive hot for titanium in aluminium, linear in ΔT, equal to an FE model with thermal initial strains.
13. **Lesson joint** in a realistic range; **main** writes the figure.

### Compare

The FE solver is the established method here, and it is your own; the comparison that matters is with mathematics it does not share. `scipy.integrate.quad` integrates the cone compliance and the tapered bar directly, without the partial fractions or the elements, and agrees to 10⁻¹⁰. The joint diagram is checked by an FE model that knows nothing about joint constants: two springs in parallel, one of them too short. That both routes give F_b = 23 700.8 N at the limit load is the evidence that (10) is right, rather than merely remembered.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **FE convergence.** The tapered bar with 2 to 32 elements. Predict the order.
2. **Joint parameter study.** The joint table of 2.7 for the titanium bolt and for a 4130 bolt of the same size. Predict which bolt loses more preload when cold, and which has the higher joint constant.
3. **Load sweep.** F_b and F_m against P from 0 to 40 kN at the three temperatures (the right panel of `bolted_joint.png`). Predict the slopes before and after separation.
4. **Torque scatter.** Draw 100 000 nut factors from N(0.20, 0.025) with `numpy.random.default_rng(507)`, compute the preloads from T = 32 N·m, and count the joints whose cold separation load is below the 18 kN ultimate load and whose hot bolt stress is above yield at limit or above ultimate strength at ultimate. Predict whether either fraction is below 1 %.

Raw data format: CSV tables `n_elem, delta_mm, rel_error` (1); `bolt, case, preload_kN, P_sep_kN, Fb_limit_kN, stress_limit_MPa, Fb_ult_kN, clamp_ult_kN` (2); `P_kN, T_case, Fb_kN, Fm_kN` (3); `seed, n, frac_sep_cold_ult, frac_yield_hot_limit, frac_ult_hot_ult` (4); plus `bolted_joint.png`.

What you should find: (1) order 2 (errors falling by a factor 3.4 at the coarsest step and 4.0 at the finest, approaching the asymptotic 4 per halving); (2) the steel bolt has C = 0.445 against 0.308 and loses 4.39 kN cold against 4.00 kN for titanium, because its stiffness in (12) outweighs its smaller mismatch; (3) slope C before separation and 1 after it for the bolt, −(1 − C) and 0 for the clamp; (4) **4.19 % separate below ultimate when cold, 3.55 % exceed Fty at limit when hot and 2.47 % exceed Ftu at ultimate when hot**. Neither fraction is acceptable for a primary joint. The deterministic window explains why: cold, the joint needs F_i ≥ (1 − C) × 18 kN + 4.00 kN = 16.45 kN; hot, it needs F_i ≤ 880 MPa × A_t − 2.67 kN − C × 12 kN = 25.85 kN. The window from 16.45 to 25.85 kN is ±22 % around its middle, narrower than the ±37 % (three standard deviations) of torque control. Halving the scatter (a standard deviation of 0.0125, which is what a more precise preload method, such as tightening to a measured bolt stretch, would have to deliver) closes the gap: 0.03 % and 0.01 %.

## 5. Validate and improve

| quantity | closed form | FE / quadrature | difference |
|---|---|---|---|
| fixed–fixed bar, u, N₁, N₂ (ΔT = 0 and 40 K) | (5) | FE | ≤ 2.7 × 10⁻¹² |
| tapered bar δ | (2): 1.26429 mm | FE, 32 elements: 1.26398 mm | 2.47 × 10⁻⁴, order 2 |
| constrained bar | −EαΔT | FE, every element | 10⁻¹² relative |
| frustum stiffness | (7) | `quad` | 10⁻¹⁰ |
| layered members | sliced frusta | `quad` over the double cone | 10⁻⁸ |
| joint at limit load, F_b, F_m | (10): 23 700.8, 11 700.8 N | parallel-spring FE with misfit 112.43 μm: 23 700.8, 11 700.8 N | 10⁻¹⁰ |
| thermal preload change | (12): −4.00 kN cold | FE with thermal initial strains | 10⁻¹⁰ |

Every closed form is confirmed by a method that does not share its derivation. The model's real uncertainty is not in this table: it is the cone angle (±5° changes k_m by about 12 % and C by 0.025), embedding and relaxation (a few per cent of preload), and above all the nut factor. **The engineering answer**: the lesson joint meets both requirements at its nominal preload, with a cold separation load of 23.1 kN against 18 kN and a hot bolt stress of 720 MPa against 880 MPa at limit; but under torque control with 12.5 % nut-factor scatter about 4 % of joints fail each requirement. Recommend preload control by bolt stretch (or by a torque-angle method), or a larger bolt that widens the window.

What to improve: a measured nut factor for the actual fastener, lubricant and finish; a cone angle calibrated by an axisymmetric FE model (Extension 3); temperature-dependent properties. What not to improve: the FE mesh, which is already exact for the bars and converged for the taper.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; derivations of (5), (7), (10)–(12); `bolted_joint.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question, stating the Monte Carlo assumptions, and recommending a preload method.

### Adversarial self-checks

- **What if the external load doubles?** At 24 kN (limit × 2) the assembly-temperature joint is still clamped (P_sep = 28.92 kN), but the cold joint (P_sep = 23.13 kN) has separated and the bolt carries the full 24 kN. Doubling the load does not double the bolt load before separation (slope C) and does after it (slope 1): the response is piecewise linear, and superposition fails across the separation point.
- **What if ΔT goes to zero, or the materials match?** (12) gives zero change for any ΔT when α_b = α_m (tested), and zero for any materials when ΔT = 0. A steel bolt in steel plates keeps its preload in the cold; that is the case for matched materials in high-temperature engine flanges.
- **Which assumption fails first in the real world?** The nut factor: it is the largest uncertainty and the reason the Monte Carlo study exists. After it, the cone model near an edge, and embedding in the first hours of service.
- **What measurement could falsify the model?** Strain-gauge a bolt (or measure its stretch ultrasonically), tighten it, and load the joint: the slope of bolt force against external load is C, directly. A measured slope of 0.5 for the lesson joint would mean the member stiffness is half the cone prediction.
- **Does the result survive a change of method?** The closed forms survive FE and quadrature to round-off. The Monte Carlo fractions with 100 000 samples have a sampling standard error of about √(0.04 × 0.96/10⁵) ≈ 0.06 percentage points; seed 15 instead of 507 gives 4.33 % and 3.46 % (Problem 4).
- **Which decision changes if the uncertainty doubles?** With a nut-factor standard deviation of 0.05, both fractions grow to about a fifth (19.35 % and 18.59 %) and no torque can satisfy both requirements; at half the assumed scatter, 0.0125, both fall to a few hundredths of a per cent. The choice between torque control and stretch control is driven entirely by this number.

## Check yourself

1. A 1 m steel rod (E = 205 GPa) of 50 mm² carries 10 kN. How much does it stretch?
2. A bar fixed between rigid walls is heated. Why does its stress not depend on its length or area?
3. In a statically indeterminate bar, which equation is added to equilibrium, and what physical statement does it express?
4. Define the joint constant C. Which way does it move if the clamped plates are made thicker (stiffer)?
5. A joint has F_i = 20 kN and C = 0.3. What are the bolt and clamp forces at P = 10 kN, and at what P does it separate?
6. A titanium bolt clamps aluminium plates. Does the preload rise or fall when the joint is cooled, and why?
7. Why does a preloaded bolt survive fatigue better than an unpreloaded one under the same alternating load?
8. Why is a ±12.5 % scatter of the nut factor a problem for a joint with both a minimum-clamp and a maximum-stress requirement?
