# P05.11  Shear Stress in Beams

> Bending stress changes along a beam, and the difference has to be passed from fibre to fibre. That passing is shear, and in a thin wall it flows like a current: in at the free edges, along the walls, out through the web.

| | |
|---|---|
| **Learning outcomes** | 1. Derive τ = VQ/(It) from the longitudinal equilibrium of part of a beam slice, and apply it to rectangles, circles and I-sections, including the jump at a web–flange junction. 2. Define shear flow q = τt, derive the general open-section shear-flow formula for sections with a product of inertia, and apply the junction and free-edge rules to branched sections. 3. Compute the shear flow at a spar cap–web joint and the fastener pitch it requires. 4. Derive and compute the shear centre of a channel, say why a load through the web twists it, and locate the shear centre of symmetric, antisymmetric and leg-junction sections by argument. 5. Implement and test a thin-walled shear-flow solver and check it against an independent route: the shear flow built directly from the change of bending stress between two stations. |
| **Prerequisites** | `P05.10` (flexure and unsymmetric bending formulas, the spar of this lesson), `P05.09` (shear force diagrams), `P05.05` (first and second moments of area), `P05.08` (shear stress and torsion), `P01.10` (line integrals along a path), `P00.14`. Skill check: you can compute the first moment Q of a rectangle about an axis outside it. |
| **Engineering question** | The light-aeroplane spar of `P05.10` carries 13 975 N of shear at its root at limit load. What shear stress does that put in its 3 mm web, how much of the shear do the caps carry, what shear flow must the cap–web joint transfer and at what pitch must fasteners be placed if the caps were fastened on, and why does a channel-section spar twist when loaded through its web? The answers: 29.1 MPa at the neutral axis, the caps carry 5.4 %, 73.4 N/mm at the joint, a pitch of 72.6 mm for two rows of fasteners assumed good for 4 kN each, and a channel 100 × 40 × 2 mm has its shear centre 14.1 mm outside the web. |
| **Assumptions and limits** | Euler–Bernoulli beam with the bending stress of `P05.10`; shear stress uniform across the width t of the cut (good for thin walls, poor across wide flanges); thin walls (t ≪ section size), wall bending stiffness about its own midline neglected; open sections only (a closed cell such as a wing box needs one more unknown, from torsion: `P05.08`); straight wall segments; shear force applied through the shear centre (no torsion) unless stated; no shear lag, no web buckling, no stress concentrations at fasteners. |
| **Analytical work** | VQ/(It) from slice equilibrium; the parabola in a rectangle, 1.5 V/A and 4V/(3A); the I-section web stress and the junction jump; the wall-element equilibrium ∂q/∂s + t ∂σ/∂z = 0 and the general open-section formula; flange and web flows in an I-beam; cap–web joint flow and fastener pitch; the channel shear centre e = 3b²t_f/(6bt_f + ht_w). |
| **Computational deliverable** | `starter/main.py`: `vq_it`, `rect_shear_stress`, `circle_max_shear`, `i_section_shear_stress`, `thin_wall_properties`, `open_section_shear_flow`, `shear_flow_resultant`, `shear_centre`, `channel_shear_centre_closed_form`, `shear_flow_from_bending_gradient`, `fastener_pitch` (`_incoming`, `channel_segments`, `i_beam_segments`, `z_segments` and `main` are given). Tests prove: the rectangle stress integrates to V with peak 1.5 V/A; the I-section web carries more than 90 % of V and the stress jumps by b_f/t_w at the junction; thin-walled properties match hand formulas; the shear flow's resultant equals (S_x, S_y), it vanishes at free edges, is continuous at junctions and linear in the load; it matches the closed forms in the I-beam flanges and web and the independent bending-gradient route on three sections; the channel's shear centre matches the closed form and leaves no moment; the I- and Z-sections have theirs at the centre. |
| **Practical activity** | **computational experiment**: thin-wall against thick-section formulas, a shear-centre sweep for channels, a polygon-convergence study of the shear centre of a slit tube, and the fastener pitch along a spar. |
| **Equipment** | None. A computer with the environment from `P00.01`. A plastic channel or an aluminium angle from a hardware shop, loaded by hand through its web, shows the twist; nothing is measured. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `shear_flow.png`, the derivations of VQ/(It), the open-section formula and the channel shear centre, the four experiment tables, and a half-page note on the spar web and the cap–web joint. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: VQ/(It) is derived from a free-body diagram of part of a slice; all public tests pass; the shear flow is checked against the bending-gradient route and its resultant; the note gives the web stress, the joint flow and the fastener pitch with the assumed fastener strength stated. |
| **References** | `@mit-ocw-3-11-mechanics-of-materials` Roylance module "Stresses in Beams" (https://ocw.mit.edu/courses/3-11-mechanics-of-materials-fall-1999/resources/mit3_11f99_bstress/): Shear stresses; Example 2; Example 3 (required-free): VQ/(Ib) from slice equilibrium, the parabolic distribution in a rectangle and the web–flange junction of I- and T-beams. `@tudelft-alderliesten-aerospace-structures-materials` Ch. 5 "Aircraft & spacecraft structures": 5.6 Torsion box; Ch. 7 "Translating loads to stresses": 7.3 Torsional loading of a fuselage structure (required-free, read online, link only): why aircraft carry shear in thin skins and webs. |
| **Extensions** | ★ Add a closed cell: make the I-beam into a box by joining the flange tips with skins, cut it at one point, compute the open-section flow, and add the constant flow q₀ that makes the twist zero (the wing-box shear centre). ★ Compute the warping (axial displacement) of the channel under torsion and say why open sections are so much more flexible in torsion than closed ones. ★ Replace straight segments by an arc integrated with `scipy.integrate.quad` and compare the slit-tube shear centre with Experiment 3. |
| **Time** | unmeasured |

## 1. The engineering problem

In `P05.10` the spar caps carried the bending moment of the wing as a couple. But the moment changes along the span: 29.7 kN·m at the root, zero at the tip. If the bottom cap carries 156 kN of tension at the root (259.8 MPa at its centroid times 600 mm²) and almost nothing near the tip, something must keep pulling on it all along its length, and that something is the web, through the joint between web and cap. The web passes force into the caps by shear. In a riveted or bolted spar that shear goes through the fasteners, and if they are spaced too far apart they shear off; in an integrally machined spar it goes through the web–cap fillet. The same shear flow carries the wing's lift force inboard: the web is what holds the wing up against gravity in a hangar and against lift in flight.

Shear also explains a behaviour that surprises everyone the first time: a channel-section beam (a C) loaded straight down through its web does not just bend, it twists. Control-surface spars, ribs and stiffeners are often channels, and a designer who loads one through the web without allowing for that twist gets an aileron that flutters or a rib that cracks at its attachment. The point through which the load must pass to avoid the twist is the **shear centre**, and for a channel it lies outside the section altogether.

This lesson derives the shear stress from equilibrium, extends it to thin-walled sections as a shear flow, computes the spar web, the cap–web joint and the channel's shear centre, and checks the shear-flow solver against an independent route that never uses the VQ/(It) formula. Success criteria, written first: the shear flow has resultant (S_x, S_y) to round-off, vanishes at every free edge and is continuous at every junction; it agrees with the bending-gradient route; the channel's shear centre matches the closed form to 10⁻⁹; the I-section web carries more than 90 % of the shear.

## 2. Understand and derive

### 2.1 Why there is shear: a free-body diagram of part of a slice

Take a slice of beam between stations z and z + dz, where the bending moments are M and M + dM, with dM/dz = −V in the convention of `P05.10` (V the shear carried by the section, the lift outboard of it for the wing). Now cut the slice again horizontally at height y₁ and keep the part above the cut, of area A′. On its two end faces act the bending stresses σ = −My/I and −(M + dM)y/I. They do not balance:

```
         stations z               z + dz
           ┌───────────────────────────┐ ← top of section
  σ(z) ──► │      part above y₁        │ ◄── σ(z + dz)        net push ∝ dM · Q / I
           └═══════════════════════════┘ ← cut at y₁: shear τ on the cut face (width t)
                 τ t dz  ──►
```

The net axial force on the part is ∫_{A′} (dM y/I) dA = dM Q/I, where **Q = ∫_{A′} y dA** is the first moment of the area beyond the cut about the neutral axis. The only other axial force on the part is the shear stress τ on the horizontal cut face, of width t and length dz. Equilibrium along z gives τ t dz = |dM| Q/I, and with |dM/dz| = V,

**τ = V Q / (I t).**   (1)

This is the "VQ/It" formula (Roylance, "Stresses in Beams", Shear stresses). It is pure equilibrium plus the flexure stress; the shear stress on the horizontal cut equals the vertical shear stress at the same point because shear stresses on perpendicular planes are equal (`P05.08`). Units: N × m³ / (m⁴ × m) = N/m².

### 2.2 Solid sections

**Rectangle** b × h: Q(y) = b(h/2 − y)(h/2 + y)/2 = (b/2)(h²/4 − y²), so

**τ = (3V / 2bh)(1 − 4y²/h²),**   (2)

a parabola, zero at the top and bottom faces (which carry no shear) and **1.5 V/A** at the neutral axis, half as much again as the average V/A. Integrating τb over the depth returns V (a test). **Solid circle**: Q at the centre is 2r³/3 and t = 2r, so τ_max = **4V/(3A)**. These factors, 1.5 and 4/3, are why a shear pin or a solid axle is checked with more than the average.

### 2.3 The spar web

For the lesson's root spar (h = 180 mm, caps 60 × 10 mm, web 3 mm, I = 9.704 × 10⁶ mm⁴, V = 13 974.5 N at limit load), (1) with t = t_w in the web gives:

- at the neutral axis, Q = 60 × 10 × 85 + 3 × 80²/2 = 60 600 mm³ and **τ = 29.09 MPa**;
- at the top of the web, Q = 51 000 mm³ (the cap alone) and τ = 24.48 MPa: the web stress is nearly uniform;
- just inside the cap, the same Q is divided by b_f = 60 mm instead of 3 mm: τ = 1.22 MPa, a jump by the width ratio 20.

The average web stress V/(h_w t_w), with h_w = 160 mm the clear web depth, is 29.11 MPa, within 0.1 % of the peak; V/(ht_w) with the overall depth gives 25.88 MPa. Integrating τ over the web, **the web carries 94.6 % of V**, the caps 5.4 %. So aircraft designers size spar webs with τ ≈ V/(h_w t_w) and do not bother with the parabola. Against a typical ultimate shear strength of 331 MPa for 7075-T6 (ASM Aerospace Specification Metals data sheet; a typical value, not an allowable) the strength margin is large; webs this thin are sized by buckling instead (Problem 6).

The vertical stress that (1) gives in the flanges is small and not physically meaningful: the flange's upper and lower faces are free, so the vertical shear must vanish there, and the real shear in a flange runs horizontally, along the flange. To see it we need the thin-walled view.

### 2.4 Shear flow in thin walls

In a thin wall, the shear stress runs along the wall's midline and is nearly uniform through the thickness t. Its product with t is the **shear flow** q = τt (N/m), the force per unit length of wall. Measure arc length s along the wall. A wall element ds × dz has σ on its faces at z and z + dz and q on its edges at s and s + ds; equilibrium along z gives

**∂q/∂s + t ∂σ/∂z = 0.**   (3)

Insert the unsymmetric bending stress of `P05.10`, equation (6), with dM/dz = −S_y and dM_h/dz = −S_x, (S_x, S_y) being the shear force on the section. Integrating (3) from a point where q is known gives the **open-section shear flow**

**q(s) = q₀ − [(S_x I_xx − S_y I_xy)/D] ∫₀ˢ t x ds − [(S_y I_yy − S_x I_xy)/D] ∫₀ˢ t y ds,  D = I_xx I_yy − I_xy².**   (4)

With I_xy = 0 and S_x = 0 it is q = q₀ − (S_y/I_xx) ∫ t y ds: the thin-walled VQ/I, with the integral playing the role of Q. Two rules close the problem, and they are Kirchhoff's current law in disguise:

- **At a free edge q = 0**: the wall's free edge carries no shear, so start every integration there.
- **At a junction, flow in equals flow out**: the axial equilibrium of a tiny piece of wall around the junction.

For an open section with straight walls, x and y are linear along each segment, so the integrals in (4) are quadratic polynomials in s and can be written in closed form. The lesson's solver takes the segments in an order where every segment starts at a free edge or at a junction whose incoming flows are already known, and checks that the flow comes back to zero at the last free edge: a global equilibrium check.

### 2.5 The I-beam, and the cap–web joint

Model the spar by its centrelines: flanges of width b = 60 mm at ±h_c/2 with h_c = 170 mm, t_f = 10 mm, web t_w = 3 mm. The thin-wall I_xx = t_w h_c³/12 + 2bt_f(h_c/2)² = 9.898 × 10⁶ mm⁴, 2 % above the exact 9.704 × 10⁶ (the flange area is counted at its centreline and the web overlaps the flanges). Under S_y = V, (4) gives:

- in each half-flange, a flow whose magnitude grows **linearly** from 0 at the tip to V t_f (b/2)(h_c/2)/I = **36 001 N/m** at the junction: it runs in from the tips toward the web in the bottom flange, and out from the web to the tips in the top flange;
- in the web, which the two half-flange flows enter at the bottom junction and leave at the top one (72 002 N/m at each junction), the flow runs upward and grows as a **parabola** to **87 303 N/m** at the neutral axis, τ = 87 303/3 mm = **29.10 MPa**, against 29.09 MPa from (1);
- the resultant of the whole distribution is exactly (0, V), and its moment about the centroid is zero.

The flow that must pass from the web into each cap is V Q_cap/I, with Q_cap the first moment of the whole cap: **q = 73 444 N/m (73.4 N/mm)** with the exact I. If the caps were fastened to the web, each fastener would carry q times the pitch p. With two rows of fasteners, each assumed good for P_allow = 4.0 kN (an illustrative value; real fastener allowables come from the fastener standard and the joint's sheet thickness), and the 1.5 factor at ultimate load,

**p = n_rows P_allow / (1.5 q) = 2 × 4000 / (1.5 × 73 444) = 72.6 mm.**   (5)

Outboard, V and q fall and the pitch can open up (Experiment 4).

### 2.6 The shear centre

Every shear-flow distribution has a resultant force and a moment. The resultant is (S_x, S_y) by construction; the moment about some point fixes the **line of action** of that resultant. The **shear centre** is the point that line passes through: a shear force applied there produces exactly the flow (4) and bends the beam without twisting it. Applied anywhere else, the force is that flow plus a torque, and an open section twists under torque very easily.

**Channel** (web height h on x = 0, flanges b pointing to +x, thicknesses t_f and t_w). Under S_y, each flange carries a force F = S_y t_f b² h/(4I) (the integral of its linear flow), away from the web in the top flange and toward it in the bottom one (the flow comes in at the bottom flange tip, runs up the web and out to the top flange tip). Those two forces form a couple F h; the web's vertical force passes through the web. The resultant S_y must therefore act at a distance e on the other side of the web such that S_y e = F h:

**e = b² h² t_f / (4I) = 3b² t_f / (6b t_f + h t_w),**   (6)

```
        e        web       flange tip
   S ●←───→│══════════════════►┐     top flange: force F away from the web
     ↑      │                         (flow linear: zero at the tip, largest at the web)
   shear    │  web: force S_y up
   centre   │  (parabolic flow)
            │◄══════════════════┘     bottom flange: force F toward the web
                                       F h (couple) = S_y e
```


using I = t_w h³/12 + 2bt_f(h/2)². For a 100 × 40 × 2 mm channel, **e = 14.118 mm** behind the web, while the centroid is 8.889 mm in front of it. A load through the centroid is 23.01 mm from the shear centre and twists the channel with 23.01 N·mm per newton of load.

Three shortcuts follow from (6)'s argument: a section with an axis of symmetry has its shear centre on that axis (the I-beam, the channel's y = 0); an antisymmetric section has it at the centre of antisymmetry (the Z-section, at the middle of the web); and a section made of straight legs meeting at one point (an angle, a T) has it at that point, because every leg's flow acts along a line through it and has no moment about it (Problem 3).

### 2.7 Where the model stops

- **Thick walls and wide flanges.** (1) assumes τ uniform across t; across a wide flange it is not, and near re-entrant corners the real stress is concentrated. Experiment 1 measures the thin-wall error for the spar.
- **Closed sections.** A wing box or a fuselage is a closed cell: (4) can start nowhere with q = 0. Cut the cell, add an unknown constant flow q₀, and fix it from the twist (Extension 1); this ties to the Bredt–Batho formula of `P05.08`.
- **Torsion of open sections.** The shear centre tells you where the load must go; if it does not go there, the open section carries the torque by warping and by tiny shear stresses through its thickness, and it is very flexible.
- **Web buckling.** Thin webs buckle in shear at stresses far below the material's shear strength. Aircraft accept this and let the buckled web carry load as a diagonal tension field, which loads the stiffeners and the caps; that design method is beyond this phase (Problem 6 shows the numbers).
- **Joints.** Fasteners carry load unevenly (the end fasteners of a row take more), holes concentrate stress, and the sheet can fail in bearing before the fastener shears.

## 3. Build the computational model

Open `starter/main.py`. `_incoming` (the junction bookkeeping), `channel_segments`, `i_beam_segments`, `z_segments` and `main` are given. A segment is `((x1, y1), (x2, y2), t)` and is directed from its first point to its second.

### Step 1: solid sections

`vq_it`, `rect_shear_stress` (2), `circle_max_shear`, and `i_section_shear_stress(V, h, b_f, t_f, t_w, y)`, which uses Q of the cap plus the web above y in the web and Q of the flange above y in the flange, dividing by t_w or b_f.

### Step 2: thin-wall properties

`thin_wall_properties(segments)` sums tL for the area and centroid, then, about the centroid, uses the exact integrals along a straight segment, ∫ty² ds = tL(y₁² + y₁y₂ + y₂²)/3 and ∫txy ds = tL(2x₁y₁ + x₁y₂ + x₂y₁ + 2x₂y₂)/6.

### Step 3: the shear flow

`open_section_shear_flow(segments, Sx, Sy, n)` computes the two coefficients of (4) once, then for each segment adds the incoming flow (`_incoming`) and subtracts the two integrals, which along a straight segment from (x₁, y₁) with unit direction (u_x, u_y) are t(x₁s + u_x s²/2) and t(y₁s + u_y s²/2) (x, y from the centroid):

```python
q = _incoming(segments, k, q_end) - kx * t * (a1 * s + ux * s**2 / 2) - ky * t * (b1 * s + uy * s**2 / 2)
```

### Step 4: resultants and the shear centre

`shear_flow_resultant(segments, flows, about)` integrates q times the unit tangent for the force and q times the moment arm (x u_y − y u_x) for the moment, by Simpson's rule, which is exact for the quadratic q of a straight wall. `shear_centre` applies a unit S_y and a unit S_x and reads off x_s = M_z/S_y and y_s = −M_z/S_x. `channel_shear_centre_closed_form` is (6).

### Step 5: the independent check

`shear_flow_from_bending_gradient(segments, Sx, Sy, dz, n)` never uses (4). It evaluates the bending stress of `P05.10` at every point of the wall at two stations whose moments differ by −S_y dz and −S_x dz, integrates the stress difference times t along the wall with the cumulative trapezoid rule, and carries flows through the junctions. It is (3) solved numerically.

### Step 6: fasteners

`fastener_pitch(q, P_allow, factor)` is (5) for one row.

Running `main()` prints

```
spar root shear V = 13974.5 N (limit); I = 9.704e+06 mm^4
VQ/It in the web: neutral axis 29.09 MPa, web top 24.48 MPa, average V/(h_w t_w) 29.11 MPa, V/(h t_w) 25.88 MPa
flange just above the junction (VQ/I b_f): 1.22 MPa
shear flow at the cap-web joint VQ_cap/I = 73444 N/m (73.4 N/mm)
if the caps were fastened to the web with two rows of fasteners rated 4.0 kN each (assumed): pitch at ultimate (1.5) 72.6 mm
fraction of V carried by the web (integral of tau t_w): 0.9464
rectangle 1.5 V/A, circle 4V/(3A): 1.50, 1.3333 times V/A

thin-wall I: I_xx = 9.898e+06 mm^4 (thick formula 9.704e+06); web q max 87303 N/m -> tau 29.10 MPa; flange-tip-to-junction q 36001 N/m each
resultant (0.000, 13974.476) N, moment about centroid 7.11e-15 N m; q at the last free edge 0.00e+00 N/m

channel 100 x 40 x 2 mm: shear centre at x = -14.118 mm, y = -3.47e-15 mm; closed form e = 14.118 mm behind the web
centroid at x = 8.889 mm: a load through the centroid is 23.01 mm from the shear centre and twists the channel by 0.0230 N m per N of shear
Z section shear centre: (9.39e-16, 1.73e-15) mm
bending-gradient route vs formula, I-beam: max |dq| = 8.29e-10 N/m (9.5e-15 relative)
```

and writes `shear_flow.png`: the vertical shear stress through the spar depth with its jump at the junctions, the thin-wall shear flow drawn along the I-section walls, and the channel with its centroid and its shear centre outside the section. (The resultant 13 974.476 N is V itself, which the first line prints rounded to 13 974.5 N: Simpson's rule is exact for the quadratic flow.)

### Verify

Run `python -m pytest curriculum/phase-05/lesson-11-shear-stress-beams`. The tests:

1. **Solid sections**: the rectangle's stress integrates to V, peaks at 1.5 V/A and vanishes at the faces; VQ/(It) built by hand at the neutral axis; the circle's 4/3 from VQ/(It) with Q = 2r³/3.
2. **I-section**: web plus flanges integrate to V; the web carries 90–100 %; the junction jump equals b_f/t_w; the web peak is within 1 % of V/(h_w t_w).
3. **Thin-wall properties**: the channel's area, centroid b²/(h + 2b) and I_xx against hand formulas; the Z-section's I_xy against its two flange contributions.
4. **Shear flow equilibrium**: for an I, a channel and a Z under three load cases, the resultant equals (S_x, S_y) to 10⁻⁹, the flow starts at zero and returns to zero at the last free edge.
5. **Continuity and linearity**: flow into each junction equals flow out; doubling the load doubles the flow.
6. **Closed forms**: the I-beam flange flow is linear and the web flow parabolic, matching the hand formulas to 10⁻¹⁰; the web flow runs upward under an upward shear.
7. **Independent route**: the bending-gradient flow equals (4) on the I, the channel and the Z.
8. **Shear centre**: the channel against (6) for three geometries, on the axis of symmetry, outside the web; zero moment about it; the flange-force couple by hand; the I and the Z at the origin.
9. **Fasteners**: pitch ∝ 1/q and 1/factor.
10. **main** writes the figure.

### Compare

The bending-gradient route agrees with (4) to 10⁻¹⁴ relative on the I-beam, because for straight walls the stress difference varies linearly along each wall and the trapezoid rule integrates a linear function exactly: the two routes are the same physics by different arithmetic, and their agreement checks the arithmetic (the junction bookkeeping, the signs, the centroid) rather than a discretisation. (When you compare the two at different sets of points, interpolation error appears instead: linear interpolation of the quadratic flow between 11 points per segment gives 147 N/m differences, which vanish when both are evaluated at the same points.) For a curved wall, the polygon approximation is the discretisation; Experiment 3 measures it.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Thin-wall against thick-section formulas.** For the spar (h = 180 mm, b = 60 mm, t_w = 3 mm) and flange thicknesses 2–20 mm, compare I and the neutral-axis shear stress from the thin-wall centreline model with the exact-section formulas. Predict which error is larger. *Finding:*

   | t_f (mm) | t_f/h | I error, thin wall | τ_NA, exact (MPa) | τ_NA, thin wall (MPa) | τ error |
   |---|---|---|---|---|---|
   | 2 | 0.011 | +1.4 % | 31.82 | 31.74 | −0.24 % |
   | 5 | 0.028 | +1.9 % | 29.67 | 29.62 | −0.15 % |
   | 10 | 0.056 | +2.0 % | 29.09 | 29.10 | +0.04 % |
   | 15 | 0.083 | +1.9 % | 29.34 | 29.42 | +0.27 % |
   | 20 | 0.111 | +1.6 % | 29.85 | 30.02 | +0.57 % |

   The second moment is off by up to 2 % (the web overlaps the flanges in the centreline model) but the web shear by less than 0.6 %, because Q and I are overestimated together. For shear flow the thin-wall model is excellent even for these thick caps.
2. **Channel shear centre.** For h = 100 mm, t = 2 mm and b/h = 0.1–1.0, compute e numerically and from (6), the centroid, and the arm from the centroid to the shear centre. Predict whether e/b rises or falls with b/h. *Finding:* e = 1.875, 5.455, 14.118, 23.478, 42.857 mm for b/h = 0.1, 0.2, 0.4, 0.6, 1.0, numerical and closed form identical to the printed digits; e/b rises from 0.19 to 0.43 (it tends to 1/2 as the flanges dominate); the arm from centroid to shear centre grows from 2.71 mm to 76.19 mm. Wide-flanged channels twist hardest.
3. **The slit tube.** A thin circular tube of radius R = 50 mm slit along its length (the shape of a spacecraft's extendable STEM boom) has its shear centre at 2R from the centre, on the side away from the slit. Approximate it by n straight segments and compute the shear centre. Predict the convergence order. *Finding:* 1.50R with 4 segments, 1.892R with 8, 1.974R with 16, 1.9936R with 32, 1.9984R with 64, 1.9996R with 128: the error falls by a factor of 4 per doubling, order 2, the error of replacing an arc by chords. The shear centre is a full radius *outside* the tube: a slit tube loaded through its centre twists violently, which is why STEM booms are loaded through end fittings designed for it.
4. **Fastener pitch along the span.** Using the elliptic shear of `P05.10` with the root section held constant, compute V, the joint flow and the pitch at y = 0–4.5 m. Predict the pitch at mid-span. *Finding:* 72.6 mm at the root, 97.2 at 1 m, 143.9 at 2 m, 255.0 at 3 m, 697.7 at 4 m, 1942 at 4.5 m; the web stress falls from 29.1 to 1.1 MPa. In practice the pitch is capped (a maximum pitch keeps the sheet from buckling between fasteners and seals the joint), so the outboard pitch is set by rules, not by shear.
5. **The Z-section's flange flow.** Under a vertical shear of 1000 N on the Z-stringer of `P05.10` (centreline 40 × 19 × 2 mm), plot the flow along the top flange. Predict whether each flange carries a net force. *Finding:* the flange flow starts at zero, reaches +3466 N/m, changes sign about two-thirds of the way to the web, and arrives at −10 401 N/m; its net force is zero (2.8 × 10⁻¹⁴ N). The web carries all 1000 N with a peak of 32 299 N/m, and the shear centre is at the middle of the web.

Raw data format: CSV tables `tf_mm, I_err, tau_exact_MPa, tau_thin_MPa` (1); `b_over_h, e_num_mm, e_closed_mm, xc_mm, arm_mm` (2); `n_segments, xs_over_R, error` (3); `y_m, V_N, q_N_per_m, pitch_mm, tau_web_MPa` (4); `s_mm, q_N_per_m` (5); plus `shear_flow.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| web τ at neutral axis | VQ/(It) = 29.09 MPa | thin-wall 29.10 MPa | 0.04 % | centreline model |
| web share of V | > 90 % | 94.6 % | — | — |
| resultant of q | (0, V) | (0, 13 974.476) N = (0, V) | < 10⁻⁹ | — |
| q at the last free edge | 0 | 0 | — | — |
| bending-gradient route | (4) | 9.5 × 10⁻¹⁵ relative | round-off | same physics |
| channel shear centre | (6): 14.118 mm | 14.118 mm | < 10⁻⁹ | — |
| slit-tube shear centre | 2R | 1.9996R at 128 segments | order 2 | chords for an arc |
| I and Z shear centre | at the centre by symmetry | 10⁻¹⁵ mm | round-off | — |

The engineering answer: **the spar web carries 94.6 % of the root shear at a nearly uniform 29 MPa; the cap–web joint must pass 73.4 N/mm, which two rows of 4 kN fasteners carry at a pitch of 72.6 mm at the root (ultimate factor included); a 100 × 40 × 2 mm channel loaded through its web is loaded 14.1 mm off its shear centre and twists.** The web stress is low against the material's shear strength; the web is sized by buckling, which the next improvement must address.

What to improve: the closed cell (Extension 1), because a real wing box is closed and its shear centre and torsional stiffness come from that; the shear buckling of the web; the uneven loading of fasteners. What not to improve: the thick-section correction to the web stress, which Experiment 1 shows is under 0.6 %.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (4) and (6); `shear_flow.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on the spar root's web stress, the cap–web joint flow and the fastener pitch, with the fastener strength you assumed and how you would verify it.

### Adversarial self-checks

- **What if the load doubles?** Everything in this lesson is linear in V: the web stress, the joint flow and the fastener load double, and the pitch halves (36 mm). The shear centre does not move: it is a property of the section, not of the load.
- **What if the web thickness goes to zero?** The web stress V/(h_w t_w) grows without bound; long before that the web buckles. The shear centre of a channel moves toward e = 3b/6 = b/2 (the limit of (6) as t_w → 0), the farthest it can be.
- **Which assumption fails first in the real world?** Web buckling for a thin web; then uneven fastener loading at the ends of a joint; then shear lag in wide caps or skins near load introductions.
- **What measurement could falsify the model?** A strain-gauge rosette on the web at the neutral axis: the shear strain should be τ/G = 29.1/26.9 GPa = 1.08 × 10⁻³ at limit load, with zero direct strain there. For the channel, a dial gauge on each flange tip under a load through the web: unequal deflections are twist; moving the load to 14 mm behind the web should make them equal.
- **Does the result survive a change of method?** The closed-form segment integrals and the bending-gradient route agree to round-off; the thin-wall and exact-section web stresses agree to 0.6 % over a tenfold range of cap thickness.
- **Which decision changes if the uncertainty doubles?** The fastener pitch is proportional to the assumed fastener strength; if that is uncertain by ±20 % instead of ±10 %, the pitch would be set from the lower bound, about 10 % closer, adding fasteners (and mass and cost) along the whole root region.

## Check yourself

1. In τ = VQ/(It), what area is Q the first moment of, and about which axis?
2. What is the peak shear stress in a solid rectangle in terms of V/A, and where is it?
3. Why does the vertical shear stress jump at the junction between a web and a flange, and by what factor?
4. Define shear flow and give its units. What two rules fix the shear flow at free edges and junctions?
5. A spar web 160 mm deep and 3 mm thick carries 14 kN of shear. Estimate the web shear stress in one line.
6. What is the shear centre, and what happens to a channel loaded through its web?
7. Where is the shear centre of an equal-legged angle, and why?
8. Fasteners rated 4 kN each, one row, joint flow 50 N/mm, factor 1.5: what pitch?
