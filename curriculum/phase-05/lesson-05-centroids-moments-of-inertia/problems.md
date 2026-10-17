# P05.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Draw each section with its axes and centroid, state units on every numerical answer (mm, mm², mm⁴ are convenient), and give the sign of every product of area.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a T-section by hand)

A T-section has a flange 60 mm wide and 8 mm thick on top of a web 6 mm thick. The overall height is h = 40 + s mm. Put the origin at the bottom of the web, x horizontal.

(a) Set up a composite table (A, ȳ, Aȳ, I_c, d, Ad²) and find the centroid and Ix.
(b) Find Iy (no table needed: both parts are centred on the y-axis) and say why Ixy = 0.
(c) Find the section moduli to the top and bottom fibres and the radius of gyration r_x. Which fibre will see the larger bending stress?
(d) Check with `composite_section` and with `polygon_section` on an 8-vertex outline.

Answer form: a table, ȳ, Ix, Iy, two section moduli, r_x.

## 2. (Derivation: the three tools of the lesson)

(a) From Green's theorem, A = ½ ∮ (x dy − y dx). Evaluate it along a straight edge from (x_i, y_i) to (x_{i+1}, y_{i+1}) and show that it gives A = ½ Σ c_i with c_i = x_i y_{i+1} − x_{i+1} y_i. Then split the polygon into triangles (origin, vertex i, vertex i + 1) and show that Q_y = ⅙ Σ (x_i + x_{i+1}) c_i.
(b) Prove the parallel-axis theorem for Ix and for Ixy, and explain why it fails if the starting axis is not centroidal.
(c) Derive the rotation formulas (4) of the lesson, show that Ix + Iy and Ix Iy − Ixy² are invariant, and find the angle θ_p at which I_uv = 0. Then show that the second moment about an axis with unit normal n is nᵀCn with C = [[Iy, Ixy], [Ixy, Ix]].

Answer form: three derivations.

## 3. (Unfamiliar scenario: a helicopter blade D-spar)

The leading edge of a helicopter rotor blade is often a D-shaped spar: a semicircular nose and a flat rear web. Model one (dimensions illustrative): outer radius 30 mm, centre at the origin, nose toward −x, flat back at x = 0; wall thickness 2.5 mm round the nose; rear web 2 mm thick (occupying −2 ≤ x ≤ 0). Inside the nose, a tungsten-alloy balance rod of diameter d = 6 + s/5 mm sits with its centre at x = −20 mm.

(a) Build the spar as a polygon with a hole (sample each arc with 2000 points) and compute A, the centroid and Ix, Iy. How do the results change with 50 and 200 points?
(b) Check A and x̄ with a thin-walled estimate: a semicircular arc of mid-line radius r_m (centroid 2r_m/π from the diameter) and a web of length 2r_m.
(c) The spar is 2024-T3 (typical density 2 780 kg/m³) and the rod a tungsten heavy alloy (typical density 17 000 kg/m³; commercial grades run from about 17 000 to 18 500 kg/m³). Find the mass per unit length of each and the chordwise position of their combined centre of mass, measured from the leading edge.
(d) Blade designers put the section's centre of mass well forward, typically at or ahead of the quarter chord, to avoid flutter. Explain why a dense rod in the nose is a more efficient way to do that than thickening the nose wall.

Answer form: a section-property table with a convergence column, a thin-walled check, two masses and a position.

## 4. (Numerical, uses your code: your own Z-stringer)

A Z-stringer has web height h = 20 + s/3 mm (outside to outside), flange width b = 12 mm and thickness t = 1.2 mm.

(a) Compute A, Ix, Iy, Ixy by `polygon_section` and by `composite_section` and by `z_section_thin`. Give the thin-walled errors.
(b) Find the principal second moments and angle with `principal_axes`, and check them with `principal_axes_eig`. Explain how the eigenvectors relate to θ_p, and confirm that `rotate_axes` gives I_uv = 0 at θ_p.
(c) The stringer will buckle about its axis of least second moment. Give r_min.
(d) Keeping h and t, find with `minimize_scalar` the flange width that maximises I_min/I_max, and report θ_p and I_min there.

Answer form: a table of three methods, principal values and angle, r_min, an optimum flange width with two numbers.

## 5. (Code: idealise a wing box)

Use `wing_box_parts` with the lesson's box (400 × 120 mm, webs 2 mm, four 300 mm² caps, four 60 mm² stringers per skin) but skins t = 1.0 + s/30 mm thick.

(a) Compute Ix for the full, booms-only and lumped models, their errors, and the mass per metre in 2024-T3.
(b) What fraction of the full Ix comes from the two skins?
(c) Vary the number of stringers per skin (0, 2, 4, 8). How do Ix, the lumped model's error and Ix per kg/m change?

Answer form: three second moments with errors, a mass, a fraction, a table of four rows.

## 6. (Judgment: how good are nominal section properties?)

An analyst computes the Z-stringer of the lesson (25/15/1.6 mm) with the thin-walled formulas and the I-spar from its nominal drawing dimensions.

(a) Compare the thin-walled error of the Z with the effect of a ±0.1 mm thickness tolerance.
(b) Two 4.8 mm rivet holes are drilled through the I-spar's top flange at one station. Estimate the change in the centroid and Ix at that station. Does the hole change the spar's overall stiffness?
(c) The extrusion has fillet radii at the web–flange junctions. Without calculating, say whether they raise or lower A and Ix, and by roughly how much compared with the other effects.
(d) Which of these would you account for in a preliminary design, and which in a stress check at a fastener row?

Answer form: two estimates and a short engineering argument.
