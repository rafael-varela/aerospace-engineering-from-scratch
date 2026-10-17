# P05.08  Torsion of Shafts

> A shaft carries power as torque times speed. Make it hollow and it barely notices; cut a slit along it and it forgets how to resist at all.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the torsion formulas τ = Tr/J and φ = TL/GJ for circular shafts from compatibility (sections rotate rigidly), Hooke's law in shear and equilibrium of the cut section, and relate torque to power and speed. 2. Size solid and hollow shafts against strength at limit and ultimate load and against a twist requirement, and compare their mass. 3. Derive the Bredt–Batho relations for a single-cell thin-walled closed section, apply them to a wing box, and compare its stiffness with the same walls cut open. 4. Implement and test an independent check: Prandtl's stress function solved by finite differences, converging at second order to the exact series for a solid rectangle and approaching Bredt–Batho for a thin square tube. |
| **Prerequisites** | `P05.05` (polar second moment J, the shoelace area of a polygon), `P05.06` (shear stress, shear strain, G = E/2(1 + ν)), `P05.07` (statically indeterminate members by compatibility), `P04.04` (rotation; power P = Tω), `P02.11` (the Poisson equation), `P03.11` and `P03.09` (finite differences and sparse solvers), `P01.12` (Green's theorem). Skill check: you can compute the polar moment of a hollow circle and solve a 5-point Laplacian on a square grid. |
| **Engineering question** | A turboshaft engine delivers 1400 kW at 20 900 rpm (a typical class of helicopter engine; the numbers are illustrative) through a 600 mm shaft to the gearbox. Taking limit torque as 1.25 × rated (transient over-torque) and ultimate as 1.5 × limit, and requiring the twist at limit torque to stay below 5° (an assumed drive-train stiffness requirement), what is the lightest Ti-6Al-4V shaft, and how much lighter is it hollow than solid? And for the wing, how much stiffer in torsion is a closed 400 × 120 mm box than the same walls cut open? The answers: rated torque 639.7 N·m, limit 799.6 N·m, ultimate 1199.4 N·m; strength alone needs a 22.05 mm solid shaft, but twist governs: 33.90 mm solid (4.00 kg/m) or 38.67/30.94 mm hollow (1.87 kg/m, 47 % of the solid mass) with a limit shear stress of 119 MPa; the closed box is 4886 times stiffer than the open one. |
| **Assumptions and limits** | Linear elasticity, small twist; circular shafts: plane sections stay plane and rotate rigidly (exact for circles, false for any other shape); uniform torque between stations; no stress concentration at shoulders, splines or keyways; thin-walled closed sections: shear stress uniform through the thickness and the wall thin compared with the cell size; Prandtl's formulation for non-circular solid and hollow sections; no warping restraint at the ends; no torsional buckling of thin tubes, no whirling (bending critical speed), no fatigue. Material values are typical, not design allowables. |
| **Analytical work** | Power, speed and torque; the kinematics of a twisted circular shaft; τ = Tr/J and φ = TL/GJ; solid and hollow J; stepped shafts and the fixed–fixed shaft by compatibility; sizing formulas for strength and twist; the shear flow in a thin closed wall and the Bredt–Batho relations q = T/2A_m and dθ/dx = T∮ds/t/(4A_m²G); the open-section constant Σbt³/3; Prandtl's stress function, ∇²φ = −2Gθ, T = 2∫φ dA, and the link between the stress function in a thin wall and the shear flow. |
| **Computational deliverable** | `starter/main.py`: `shaft_torque`, `polar_moment`, `torsion_shear_stress`, `twist_angle`, `stepped_shaft_twist`, `fixed_fixed_shaft`, `size_solid_shaft`, `size_hollow_shaft`, `size_for_twist`, `shaft_mass_per_length`, `polygon_area`, `bredt_batho`, `thin_tube_J`, `open_section_J`, `rectangle_torsion_constant`, `prandtl_fd` (`MATERIALS`, `shear_allowables`, the design case and `main` are given). Tests prove: P = Tω; J by quadrature; the shear stress integrates back to T (equilibrium of the cut section); surface shear strain rφ/L = τ/G; the indeterminate shaft obeys equilibrium and zero total twist; sized shafts sit exactly at the allowable; hollow tends to solid; the Bredt–Batho shear flow is constant and its moment about any point is T; the result does not depend on the starting vertex or direction; a 4000-gon reproduces the thin circular tube, whose ratio to the exact J is 1/(1 + (t/2r_m)²); the slit tube is 3(r/t)² softer; Prandtl FD converges at order 2.0 ± 0.1 to the series, is zero on the boundary and symmetric, and approaches Bredt–Batho as the wall thins. |
| **Practical activity** | **computational experiment**: a shaft-sizing trade (material and diameter ratio against mass, stress and twist), a closed-against-open section comparison, and a method comparison between Bredt–Batho and a finite-difference Prandtl solution as the wall thickens. |
| **Equipment** | None. A computer with the environment from `P00.01`. A cardboard tube, slit along its length and twisted next to an intact one, shows the open-section effect in the hand; nothing is measured with it. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `torsion.png`, the derivations of τ = Tr/J, the Bredt–Batho relations and the FD discretisation, the experiment tables, and a half-page note recommending a shaft with its governing criterion and stating the closed-against-open stiffness ratio. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: τ = Tr/J and Bredt–Batho are derived from compatibility, Hooke's law and equilibrium (not quoted); all public tests pass; the FD convergence order is reported; the shaft note names the governing criterion and the assumptions behind the limit factor and the twist requirement. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Shear and Torsion": Shearing stresses and strains; Statics - Twisting Moments; Torsional Stresses and Displacements; Noncircular sections: the Prandtl membrane analogy (required-free, https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_torsion/): the derivation of τ = Tr/J and θ = TL/GJ, and Example 4, which sizes a shaft from power and speed. `@mit-ocw-16-001-unified-materials-structures` lecture notes "Torsion" (required-free). `@libretexts-baker-haynes-engineering-statics` 10.5 Polar Moment of Inertia (required-free): J for solid and hollow sections. `@tudelft-alderliesten-aerospace-structures-materials` 7.3 Torsional loading of a fuselage structure (required-free, link only: https://books.open.tudelft.nl/home/catalog/book/197): closed thin-walled sections in an airframe. |
| **Extensions** | ★ Add a check against torsional buckling of the thin hollow shaft and the first bending critical speed (whirl) of the shaft between its bearings, and see whether d_i/d_o = 0.9 survives them. ★ Extend `bredt_batho` to a two-cell wing box (two unknown shear flows, equal twist of both cells) and compare with `prandtl_fd` on a rectangle with two holes. ★ Use `prandtl_fd` to find the maximum shear stress (the steepest slope of φ) of a solid rectangle and compare with the series value at the middle of the long side. |
| **Time** | unmeasured |

## 1. The engineering problem

A helicopter turboshaft engine turns its power turbine at about 20 000 rpm and hands the power to the main gearbox through a short output shaft. The shaft carries the full engine power as torque, P = Tω: 1400 kW at 20 900 rpm is 639.7 N·m. That is modest (a car's driveshaft in first gear sees more), but the shaft must be light, it must survive transient over-torque, and it must not be too soft: the engine governor, the gearbox and the rotor form a torsional spring–mass system, and a soft shaft lowers its natural frequency into the range where the governor can excite it. Strength and stiffness are separate requirements, and either can govern.

The same physics, in a different shape, decides the torsional stiffness of a wing. A wing twisted by aerodynamic and aileron loads resists with the closed box formed by its skins and spars. Close the box and it is very stiff; cut a slit in it (an access panel, a door, a slot for a control surface) and the stiffness collapses, by a factor of thousands. Aileron reversal and flutter, the aeroelastic failures of `P13`, start from this number.

This lesson builds both. Part 2 derives the circular-shaft formulas from compatibility, Hooke's law in shear and equilibrium, then Bredt–Batho for a thin closed section, then Prandtl's stress function, which covers any section and gives the independent check. Part 3 implements them. Part 4 sizes the shaft and compares closed with open sections and Bredt–Batho with a finite-difference solution. Success criteria, written first: the stress distribution integrates back to the applied torque; every sized shaft sits exactly at its allowable; the FD solution converges at second order to the exact series on a rectangle; and Bredt–Batho agrees with the FD solution to within a few per cent for a wall 5 % of the section width.

## 2. Understand and derive

### 2.1 Power, speed and torque

A torque T turning at angular speed ω does work at the rate **P = Tω** (`P04.04`). With the speed in rev/min, ω = 2πN/60. For the turboshaft, ω = 2188.6 rad/s and T = 1.4 × 10⁶ W / 2188.6 rad/s = 639.7 N·m. Doubling the speed at the same power halves the torque, which is why the power turbine drives a fast, light shaft and the gearbox, not the turbine, makes the large rotor torque.

### 2.2 Circular shafts: compatibility, Hooke and equilibrium

**Compatibility.** Twist a circular shaft and, by symmetry, every cross-section stays plane and circular and rotates rigidly about the axis: nothing distinguishes one diameter from another, so nothing can warp. If two sections dx apart rotate by dφ relative to each other, a line on the surface at radius r tilts by the shear strain

γ = r dφ/dx.

The shear strain grows linearly from zero on the axis to its maximum on the surface.

**Constitutive law.** Hooke's law in shear (`P05.06`): τ = Gγ = G r dφ/dx.

**Equilibrium.** The shear stresses on the cut face must add up to the internal torque: T = ∫ τ r dA = G (dφ/dx) ∫ r² dA = G J dφ/dx, with the **polar second moment** J = ∫r² dA (`P05.05`). Eliminating dφ/dx:

**τ = T r / J,   φ = T L / (G J)**   (1)

for a uniform shaft of length L. For a solid circle J = πd⁴/32; for a hollow one **J = π(d_o⁴ − d_i⁴)/32**. The maximum stress is at the outer surface, τ_max = 16T/(πd³) for a solid shaft. Material near the axis carries little stress and contributes little to J (it enters as r³ dr), which is why hollow shafts are efficient: a tube with d_i/d_o = 0.8 removes 64 % of the area but only 41 % of J.

A stepped shaft is segments in series: φ = Σ TᵢLᵢ/(GᵢJᵢ). A shaft fixed at both ends with a torque applied between them is statically indeterminate, solved by compatibility exactly as the fixed–fixed bar of `P05.07`: the total twist is zero, so the torque divides in proportion to the segments' torsional stiffnesses GJ/L. A Ti-6Al-4V shaft of a 200 mm solid 24 mm segment and a 400 mm hollow 30/24 mm segment, loaded with the 799.6 N·m limit torque at the junction, sends T₁ = 464.7 N·m into the solid segment and 334.9 N·m into the hollow one, and the junction rotates 3.856°.

### 2.3 Sizing a shaft

Setting τ_max = τ_allow in (1) gives the smallest diameters:

**d = (16T/(π τ_allow))^(1/3),   d_o = (16T/(π τ_allow (1 − k⁴)))^(1/3),  k = d_i/d_o,**   (2)

and the twist requirement φ ≤ φ_allow gives **d_o = (32TL/(π G φ_allow (1 − k⁴)))^(1/4)**. The design takes the largest diameter of all the criteria.

The lesson's criteria: at limit torque (1.25 × rated = 799.6 N·m) the surface shear stress may not exceed the shear yield strength, estimated as Fsy = Fty/√3 (the von Mises estimate derived in `P05.14`); at ultimate torque (1.5 × limit = 1199.4 N·m) it may not exceed the ultimate shear strength, estimated as Fsu = 0.6 Ftu (a common approximation for ductile metals); and the twist at limit over 600 mm may not exceed 5°. For typical Ti-6Al-4V (Fty = 880 MPa, Ftu = 950 MPa, G = E/2(1 + ν) = 42.4 GPa) those give Fsy = 508 MPa and Fsu = 570 MPa. The ultimate criterion needs a 22.05 mm solid shaft, the yield criterion 20.01 mm, and the twist criterion **33.90 mm**: stiffness governs. A strength-sized 22 mm titanium shaft would twist about 28° at limit torque over 600 mm (φ scales as d⁻⁴: (33.90/22.05)⁴ × 5°), a soft torsion spring, not a shaft.

### 2.4 Thin-walled closed sections: Bredt–Batho

Now take a tube of any cross-section shape whose wall thickness t (which may vary around the perimeter) is small compared with the cell size. Assume the shear stress is uniform through the thickness and directed along the wall. Define the **shear flow** q = τt, a force per unit length of wall.

**Equilibrium of a wall element.** Cut a small element of wall between two stations along the perimeter, and two cross-sections dx apart. The only forces along the tube axis are the shear flows on the two perimeter cuts (q₁ dx and q₂ dx, by the complementary shear stresses), so q₁ = q₂: **the shear flow is the same all the way round the cell**, even where the thickness changes. Thin walls carry high stress, thick walls low.

**Equilibrium of the section.** The wall element ds carries the force q ds along the wall. Its moment about any point O is q ds × p, with p the perpendicular distance from O to the wall's line, and p ds is twice the area of the triangle from O to the element. Summed around the cell, the triangles fill the enclosed area A_m (measured to the wall's mid-line), so

**T = 2 A_m q,   q = T/(2A_m),   τ = q/t.**   (3)

**Compatibility and energy.** The twist rate follows from equating the work of the torque to the strain energy in the walls. Per unit length, ½ T dθ/dx = ∮ τ²/(2G) t ds = (q²/2G) ∮ ds/t, and with q = T/2A_m:

**dθ/dx = T ∮ds/t / (4 A_m² G),   J = 4A_m² / ∮ ds/t.**   (4)

For a thin circular tube of mid-line radius r_m, (4) gives J = 2πr_m³t; the exact π(r_o⁴ − r_i⁴)/2 is larger by the factor 1 + (t/2r_m)²: 0.25 % at t/r_m = 0.1, and 6.25 % at t/r_m = 0.5, where Bredt–Batho is 5.9 % low.

**A wing box.** A single-cell 2024-T3 box with mid-line 400 × 120 mm, skins 1.6 mm and spars 3.0 mm, carrying T = 10 kN·m: A_m = 48 000 mm², **q = 104.2 N/mm**, τ = 65.1 MPa in the skins and 34.7 MPa in the spars; ∮ds/t = 2 × 400/1.6 + 2 × 120/3.0 = 580, J = 1.589 × 10⁷ mm⁴, and the twist is 1.312° per metre.

### 2.5 Open sections

Cut the same walls open along their length and the shear flow can no longer circulate around the cell: the torque is now carried by stresses that reverse across the thickness of each wall, with a tiny lever arm of order t. For a thin strip of width b and thickness t, J = bt³/3 (the limit of the rectangle in 2.6), and for a section made of strips, **J = Σ bᵢtᵢ³/3**. The cut wing box has J = 3252 mm⁴, **4886 times less** than the closed box. A slit circular tube is softer than the closed one by **3(r_m/t)²**: 300 for r/t = 10. At the same 10 kN·m the open box would need a shear stress of 9224 MPa: it would twist until it failed long before carrying the load. This is why access holes in torsion boxes are reinforced with doublers and why a wing with a hinged leading edge carries its torque in a separate closed box.

### 2.6 Prandtl's stress function: any section, and an independent check

For a non-circular section, cross-sections warp and (1) is wrong. Prandtl's formulation handles any shape. Write the shear stresses on the section as derivatives of a **stress function** φ(x, y),

τ_xz = ∂φ/∂y,   τ_yz = −∂φ/∂x,

which satisfies equilibrium automatically. Compatibility of the warping displacement (its mixed derivatives must agree) then requires

**∇²φ = −2Gθ**  in the section,   φ = 0 on the outer boundary,   (5)

with θ the twist per unit length, and equilibrium of the section gives **T = 2∫φ dA**, so J = 2∫φ dA/(Gθ). On the boundary of a hole, φ is a constant φ_h (the hole is free of traction), and T = 2∫φ dA + 2φ_h A_hole. The stress is the slope of the φ surface; Roylance's module describes the **membrane analogy**: φ has the shape of a soap film stretched over the section and inflated by a pressure, so the stress is highest where the film is steepest, at the middle of the long sides.

Two exact results check the method. For a solid a × b rectangle (a ≥ b) the series solution of (5) gives

**J = (ab³/3)[1 − (192/π⁵)(b/a) Σ_{n odd} tanh(nπa/2b)/n⁵]**:   (6)

0.1406 a⁴ for a square, and ab³/3 (1 − 0.630 b/a) for a thin strip, so a thin strip tends to bt³/3. For a thin wall, φ rises almost linearly from 0 at the outer surface to φ_h at the inner one, so its slope across the wall is φ_h/t = τ: **the hole value φ_h is the shear flow q (per unit Gθ)**, and T = 2φ_h A_hole + (a small wall term) becomes Bredt–Batho's T = 2A_m q. Bredt–Batho is the thin-wall limit of (5).

**Finite differences.** On a square grid of spacing h, (5) is the 5-point Poisson problem of `P03.11`. The hole is handled by making every grid node in the hole and on its boundary one unknown φ_h. The cleanest way to get the right equation for it is to minimise the discrete energy

Σ_edges (φᵢ − φⱼ)²/2 − 2h² Σ_nodes φᵢ

over all unknowns (the stationarity condition of (5)'s energy functional). Its derivative with respect to an interior node gives the 5-point equation; with respect to φ_h it gives the discrete statement that the circulation of the stress round the hole equals 2Gθ times the enclosed area, the condition that makes the warping single-valued. Edges inside the hole join a node to itself and carry no energy. J is 2h²Σφ over all nodes, hole nodes included.

### 2.7 Where the models stop

- Equation (1) is exact only for circles. Shoulders, splines and keyways raise the local stress by a stress-concentration factor, often well above 1, which fatigue design must include.
- Thin hollow shafts buckle in torsion (a wrinkling of the wall) at a stress that falls as the wall thins; at d_i/d_o = 0.9 the wall is 2.2 mm on a 44 mm shaft and the check is not optional (Extension 1).
- Bredt–Batho assumes thin walls and free warping; at a rib or a fuselage frame warping is restrained, and the stresses near it differ.
- The fatigue of a turbine shaft (from bending as it spins, and torque ripple) usually governs its detail design; this lesson's static criteria are the first pass.
- The limit factor 1.25 and the 5° twist requirement are assumptions of the lesson's design case, not values from a specification; change them and the answer changes.

## 3. Build the computational model

Open `starter/main.py`. `MATERIALS` (typical values, with G from E and ν by isotropy; status as in `P05.06`), `shear_allowables`, the design-case constants and `main` are given.

### Step 1: circular shafts

`shaft_torque(power_W, rpm)` returns T in N·mm (the lesson works in mm, N and MPa). `polar_moment`, `torsion_shear_stress`, `twist_angle` and `stepped_shaft_twist` are (1). `fixed_fixed_shaft` is the compatibility solution of 2.2: φ = T₀/(k₁ + k₂), T₁ = k₁φ, T₂ = −k₂φ.

### Step 2: sizing

`size_solid_shaft`, `size_hollow_shaft` and `size_for_twist` invert (1) as in (2); `shaft_mass_per_length` is ρ × area.

### Step 3: thin-walled sections

`polygon_area` is the shoelace formula of `P05.05`. `bredt_batho(T, vertices, t, G)` takes the mid-line corners and a thickness per wall, and returns A_m, q, τ per wall, the twist rate and J from (3) and (4). `thin_tube_J` is 2πr_m³t; `open_section_J` is Σbt³/3.

### Step 4: Prandtl by finite differences

`rectangle_torsion_constant` is (6) summed to 50 odd terms. `prandtl_fd(a, b, h, hole)` builds the energy's matrix from every horizontal and vertical grid edge with `scipy.sparse`:

```python
dof[free] = np.arange(n_free)              # interior nodes of the material
dof[in_hole] = n_free                      # one unknown for the hole and its boundary
for each edge (p, q):                      # p, q are dof indices, -1 on the outer boundary
    K[p, p] += 1; K[q, q] += 1; K[p, q] -= 1; K[q, p] -= 1   # only for p, q >= 0 and p != q
f = 2 h^2 per interior node; 2 h^2 x (number of hole nodes) for the hole unknown
phi = spsolve(K, f);  J = 2 h^2 sum(phi over all nodes)
```

(edges with one end on the outer boundary add only the diagonal term; edges inside the hole are skipped).

Running `main()` prints the sizing table, of which the titanium rows are

```
material    ratio  Fsy   Fsu   d_o (limit/yield)  d_o (ult)  d_o (twist 5 deg)  chosen d_o  d_i   kg/m   tau_lim  twist (deg)
Ti-6Al-4V   0.0   508   570      20.01            22.05        33.90            33.90   0.00  3.999     105   5.000
Ti-6Al-4V   0.6   508   570      20.96            23.09        35.10            35.10  21.06  2.743     108   5.000
Ti-6Al-4V   0.8   508   570      23.86            26.28        38.67            38.67  30.94  1.873     119   5.000
Ti-6Al-4V   0.9   508   570      28.56            31.47        44.27            44.27  39.84  1.295     136   5.000
```

then the 4130 rows, the indeterminate shaft, the thin-tube ratios, the wing box and its open counterpart, the FD convergence table and the square-tube comparison, and writes `torsion.png` (mass against d_i/d_o for both materials; τ against r for the solid and hollow shafts; the stress function of a square tube).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-08-torsion-of-shafts`. The tests:

1. **Power**: T ω = P; doubling the speed halves the torque.
2. **Equilibrium of the section**: J by `quad` of 2πr³; ∫τ r dA = T to 10⁻¹²; τ linear in r and zero on the axis.
3. **Kinematics**: twist linear in T; stepped twists add; the surface strain rφ/L equals τ/G.
4. **Indeterminate shaft**: T₁ − T₂ = T₀, zero total twist, and a rigid segment taking everything.
5. **Sizing**: every sized shaft has τ_max = τ_allow (or φ = φ_allow); hollow at k = 0 is solid; hollow is lighter.
6. **Allowables**: Fty/√3 and 0.6 Ftu; G consistent with an isotropic ν.
7. **Bredt–Batho**: q = τt on every wall; the moment of the shear flow about an arbitrary point equals T; A_m and the twist rate.
8. **Invariance**: rotating the vertex list or reversing it changes nothing.
9. **Thin tube**: a 4000-gon reproduces 2πr_m³t; the ratio to the exact J is 1/(1 + (t/2r_m)²) and rises to 1 as t → 0.
10. **Open section**: the slit tube is 3(r/t)² softer; the wing box more than 1000 times.
11. **Series (6)**: 0.1406 for a square, the thin-strip limit, symmetry in a and b.
12. **Prandtl FD**: order 1.9–2.1 against (6); φ ≥ 0, zero on the boundary, symmetric.
13. **Hollow FD**: φ flat across the hole and maximal there; Bredt/FD rising toward 1 as the wall thins, within 4 % at t = 0.05 a; a tube is softer than the solid square.
14. **main** writes the figure.

### Compare

The FD solver and Bredt–Batho share nothing but the physics: one solves a partial differential equation on a grid, the other applies equilibrium to a thin wall. The series (6) is a third, independent route for the solid rectangle. `scipy.sparse.linalg.spsolve` is the established sparse direct solver of `P03.08`/`P03.09`; the conjugate-gradient method of `P03.09` would also work, since the energy matrix is symmetric positive definite.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Shaft trade.** For Ti-6Al-4V and 4130, and d_i/d_o = 0, 0.6, 0.8, 0.9, compute the diameters required by each criterion, the governing one, the mass per metre and the limit stress. Predict which criterion governs and which material is lighter.
2. **Closed against open.** The wing box of 2.4 closed and cut open: J, the twist rate and the stress at 10 kN·m. Predict the ratio from 3(r/t)² with an equivalent radius.
3. **Bredt–Batho against Prandtl FD.** Square tubes of outer side 1 with walls t = 0.05, 0.1 and 0.2, solved with h = t/8; and the FD convergence on a solid 2 × 1 rectangle against (6). Predict how the Bredt/FD ratio depends on t, and the FD order.

Raw data format: CSV tables `material, ratio, d_yield, d_ult, d_twist, d_chosen, kg_per_m, tau_limit_MPa` (1); `section, J_mm4, twist_deg_per_m, tau_MPa` (2); `t, J_fd, J_bredt, ratio` and `h, J_fd, error` (3); plus `torsion.png`.

What you should find. (1) Twist governs every design; the closest call is 4130 at k = 0.9 (37.84 mm for twist against 36.13 mm for yield, 5 % apart); titanium is lighter at every ratio (1.873 against 2.425 kg/m at k = 0.8) despite its lower modulus, because its density is 56 % of steel's; the hollow titanium shaft at k = 0.8 weighs 0.469 of the solid one. The limit stresses (105–136 MPa in titanium) are far below the 508 MPa yield estimate: once stiffness governs, strength margins are large, which helps fatigue. (2) J falls from 1.589 × 10⁷ to 3252 mm⁴ (4886 times) and the stress at 10 kN·m from 65.1 to 9224 MPa. (3) Bredt/FD = 0.9765, 0.9476 and 0.8704 for t = 0.05, 0.1 and 0.2: Bredt–Batho underestimates the stiffness, by roughly half the ratio t/a, because the material in the corners and the thick-wall effect add stiffness the mid-line formula leaves out; FD errors on the rectangle fall by 4 per halving (orders 1.99, 2.00, 2.00).

## 5. Validate and improve

| quantity | closed form | independent route | difference |
|---|---|---|---|
| section equilibrium ∫τr dA | T | `quad` | 10⁻¹² |
| thin circular tube J | 2πr_m³t | exact π(r_o⁴ − r_i⁴)/2 | factor 1/(1 + (t/2r_m)²): 0.99751 at t/r_m = 0.1 |
| solid 2 × 1 rectangle J | series (6): 0.457363 | Prandtl FD, h = 0.0125: 0.457232 | 1.32 × 10⁻⁴, order 2 |
| square tube, t = 0.05 | Bredt–Batho 0.04287 | Prandtl FD 0.04390 | 2.3 % |
| wing box shear flow | q = T/2A_m = 104.2 N/mm | moment of q about any point = T (test) | 10⁻¹² |
| open slit tube | J_closed/J_open = 3(r/t)² | from `bredt_batho` and `open_section_J` | exact |

Each formula is confirmed by a route that does not share its derivation, and each departure (thin-tube ratio, Bredt against FD) is the size the theory says it should be. **The engineering answer**: the turboshaft output shaft is stiffness-governed; a Ti-6Al-4V tube of 38.67 mm outside and 30.94 mm inside diameter (wall 3.9 mm) weighs 1.873 kg/m, 47 % of the 33.90 mm solid shaft that meets the same twist requirement, and runs at 119 MPa at limit torque. Going to d_i/d_o = 0.9 saves another 31 % (1.295 kg/m) but leaves a 2.2 mm wall that must pass a torsional-buckling check first. For the wing, closing the box makes it 4886 times stiffer in torsion than the same walls open.

What to improve: the buckling and whirl checks of Extension 1, stress concentrations at the splines, and the real limit-torque factor and stiffness requirement from the drive-train dynamics. What not to improve: the FD grid beyond h = t/8 for the comparison here; the FD error is already well below the Bredt–Batho modelling error it is measuring.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; derivations of (1), (3), (4) and the FD discretisation; `torsion.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a shaft with its governing criterion, its mass and its stress, and stating the closed-against-open stiffness ratio of the wing box.

### Adversarial self-checks

- **What if the torque doubles?** Stress and twist double (linear elasticity); the required diameter grows by 2^(1/3) = 1.26 for strength and 2^(1/4) = 1.19 for twist, so the governing criterion can change: the stiffness criterion grows more slowly, and at high enough torque strength would govern.
- **What if the speed goes to zero at the same power?** T = P/ω grows without bound: a shaft cannot transmit finite power at zero speed. At low speed the gearbox output shaft, not the turbine shaft, is the heavy one.
- **Which assumption fails first in the real world?** Uniform stress at the shaft's ends: splines and shoulders concentrate it. For the thin hollow shaft, torsional buckling before any of the static criteria. For the wing box, restrained warping at the root rib.
- **What measurement could falsify the model?** Twist a tube in a fixture and measure the angle with a dial gauge on an arm: the slope of torque against twist is GJ/L. A slit tube that is only 10 times softer, instead of 3(r/t)² times, means the slit is not free (it is held closed by friction or a fastener).
- **Does the result survive a change of mesh or method?** The FD J on the rectangle changes by 1.3 × 10⁻⁴ between the two finest grids; Bredt–Batho and FD differ by 2.3 % on the thinnest tube, the modelling difference, not a discretisation one. The sizing does not depend on any discretisation.
- **Which decision changes if the uncertainty doubles?** If G were uncertain by 10 % instead of 5 %, the twist-governed diameter (∝ G^(−1/4)) would move by about 2.5 % instead of 1.3 %; the design would still be stiffness-governed. If the twist requirement were relaxed from 5° to 30°, the twist diameters would shrink by 6^(1/4) = 1.57 and strength would govern: the solid shaft would fall to 22.05 mm and the k = 0.8 tube to 26.28 mm (both set by the ultimate criterion). The requirement, not the material data, is the uncertainty that matters most.

## Check yourself

1. A shaft carries 100 kW at 6000 rpm. What torque does it transmit?
2. Why do plane sections of a circular shaft stay plane, and why does that argument fail for a square shaft?
3. Write τ_max for a solid shaft of diameter d under torque T. By what factor does it change if d doubles?
4. A tube has d_i/d_o = 0.8. What fraction of the solid shaft's area does it keep, and what fraction of its J?
5. State the Bredt–Batho relation between torque and shear flow. Where in a box of varying wall thickness is the shear stress highest?
6. Why is a slit tube so much softer in torsion than a closed one? Give the ratio for r/t = 20.
7. In Prandtl's formulation, what boundary condition does the stress function satisfy on the outer boundary and on a hole?
8. A titanium shaft and a steel shaft are designed to the same twist requirement. Which is lighter, and which material property decides it?
