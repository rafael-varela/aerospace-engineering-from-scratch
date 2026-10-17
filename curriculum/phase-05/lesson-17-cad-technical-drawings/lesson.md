# P05.17  CAD & Technical Drawings I

> A drawing is a contract: every feature located once, from a datum, with a tolerance. If the shop has to guess, or can measure the same thing two ways, the drawing is wrong.

| | |
|---|---|
| **Learning outcomes** | 1. Model a flat aerospace bracket in free CAD software (FreeCAD, or any free tool) from a fully constrained sketch, pad it to thickness and add holes, keeping its dimensions as named parameters. 2. Produce a dimensioned multiview drawing with a title block, in third-angle projection, that dimensions every feature once from a datum and states units, material and general tolerances. 3. Check the drawing independently in Python: compute the part's area, mass, centre of mass and inertia tensor from its parameters by Green's theorem and prism formulas, and compare with what the CAD tool reports. 4. Decide whether a dimension scheme fully constrains the part without redundancy by the rank of its constraint Jacobian, and rebuild the geometry from the dimensions alone. 5. Run a tolerance stack-up by worst case, root-sum-square and seeded Monte Carlo, and apply the floating-fastener rule to the bolt holes. |
| **Prerequisites** | `P05.05` (centroids, second moments, parallel-axis theorem, principal axes), `P01.10` (double integrals, Green's theorem), `P01.13` (rank of a matrix, linear systems), `P03.10` (Newton and least-squares solution of nonlinear systems), `P00.14` (uncertainty, normal and uniform distributions), `P04.04` (the inertia tensor). Skill check: you can find the centroid of an L-shape by splitting it into two rectangles. |
| **Engineering question** | An actuator mounting bracket is a 4 mm 2024-T3 plate, 60 × 40 mm with one corner cut off by a sloping edge, two clearance holes for M6 bolts and a lightening hole. Can it be drawn so that a machine shop could make it from the drawing alone, with nothing to guess and nothing over-specified; what does it weigh and where is its centre of mass; and will its bolt holes line up with the mating part? The answers: 19 dimensions and geometric constraints for 19 unknowns, rank 19, no degree of freedom left and none redundant; 19.11 g with the centre of mass at (25.62, 16.94, 2.00) mm from the datum corner; and ±0.1 mm hole positions sit well inside the ±0.30 mm the floating-fastener rule allows for 6.6 mm holes and 6.0 mm bolts. |
| **Assumptions and limits** | A flat part of uniform thickness and density (a prism with through holes), which is what the mass-property formulas assume; sharp corners (no fillets; Problem 3 adds a slot made of arcs); typical 2024-T3 density, 2780 kg/m³. The constraint check is local: full rank at the nominal geometry means the dimensions fix the shape near it, not that no other, distant shape satisfies them. Tolerance stacks are linear, with independent dimensions; RSS assumes each tolerance is a ±3σ band of a normal distribution, which a shop's real process may not follow. The CAD instructions name menus as in FreeCAD 1.0; other versions and other tools differ in names but not in steps. |
| **Analytical work** | The shoelace formulas for area, first and second moments from Green's theorem; the exact circle formulas and the parallel-axis theorem for the holes; the inertia tensor of an extruded plate, including the thickness correction to the perpendicular-axis theorem; counting degrees of freedom and constraints of a sketch; the worst-case and RSS stacks and their statistical meaning; the floating-fastener rule. |
| **Computational deliverable** | `starter/main.py`: `polygon_properties`, `circle_properties`, `section_properties`, `mass_properties`, `raster_properties`, `principal_axes`, `worst_case_stack`, `rss_stack`, `monte_carlo_stack`, `floating_fastener_position_tolerance`, `constraint_residuals`, `constraint_jacobian`, `constraint_status`, `solve_sketch` (the `Bracket` dataclass, `regular_polygon`, `point_in_polygon`, `bracket_sketch`, `edge_distance_ratio` and `main` are given). Tests prove: the polygon formulas reproduce rectangles and triangles and ignore vertex order; n-gon holes converge to the exact circle at second order; the bracket's properties agree with a raster integration to 10⁻³; the inertia tensor is symmetric, positive definite, obeys J_zz = J_xx + J_yy − mt²/6 and the box formula; worst case ≥ RSS; a seeded Monte Carlo gives 3σ = RSS to 1 % and 99.73 % coverage; the lesson's dimension scheme has 0 degrees of freedom and 0 redundancies, an added dimension is reported redundant and a removed one leaves a degree of freedom; the sketch is rebuilt from the dimensions to 10⁻⁸ mm. Nothing imports a CAD package. |
| **Practical activity** | **design-build**: model, dimension and draw the bracket in FreeCAD (or another free CAD tool), export the drawing as PDF, read the CAD tool's mass properties, and check them, the dimension scheme and the tolerance stack with the tested Python module. |
| **Equipment** | A computer that runs FreeCAD 1.0 or later (free, open source, Windows, macOS, Linux; from freecad.org). Any other free CAD tool that can sketch with constraints, extrude, make a 2-D drawing and report volume and centre of mass will do. No printer needed. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | The CAD file; the drawing as PDF (three views, all dimensions, title block); a table comparing the CAD tool's volume, centre of mass and inertia with `mass_properties`; the constraint-status output for your dimension scheme; the tolerance stack table; completed `main.py` with passing tests and `bracket_check.png`. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the drawing has every feature dimensioned once from a datum, units, material and a general tolerance; the CAD and Python mass properties agree to four significant figures (or the difference is explained); the constraint check reports 0 degrees of freedom and 0 redundancies for the drawing's scheme; the stack-up states which method was used and why. |
| **References** | `@freecad-manual-van-havre` "Discovering FreeCAD": The FreeCAD interface; Navigating in the 3D view; Parametric objects. "Working with FreeCAD": Modeling for product design; Generating 2D drawings; Using spreadsheets (required-free): the workflow of sketch, constrain, pad, pocket, drawing page. Its first edition targets an older FreeCAD and adds drawing dimensions through an add-on; the menu steps below are for FreeCAD 1.0. `@ford-engineering-graphics-and-design` 1. Orthographic Projection; 2. Multiview Layout and Orientation; 7. Dimensioning Standards; 8. Conventional Tolerancing; 9. Standard Fits (required-free): projection, view layout, dimensioning rules and tolerances; 4. Section Views (supplementary-free). |
| **Extensions** | ★ Drive every dimension from a FreeCAD spreadsheet and make a 3 mm and a 5 mm version from the same model; check both with `mass_properties`. ★ Add 2 mm fillets to the outer corners in the CAD model and extend `polygon_properties` with circular-arc edges (Green's theorem on an arc has a closed form); how much mass do the fillets remove? ★ Re-dimension the holes with geometric tolerancing (a position tolerance of ⌀0.3 mm to datums A, B, C) and explain why it gives the shop more room than ±0.1 mm coordinates for the same function. |
| **Time** | unmeasured |

## 1. The engineering problem

Every metal part on an aircraft starts as a drawing, or today as a CAD model with a drawing derived from it, and every one of them is checked before anyone cuts metal. The drawing must say, without ambiguity, what shape to make, how accurately, and from what material. A dimension missing means the machinist guesses; a dimension given twice (for example the overall width *and* every segment that adds up to it) means two measurements of the same thing that cannot both be exactly right, and an inspector who rejects a good part or accepts a bad one. A mass or centre-of-mass figure from the CAD tool goes into the aircraft's weight and balance; if the model has a wrong thickness or a missing hole, that error flows into every downstream analysis.

This lesson's part is an **actuator mounting bracket**: a 4 mm plate of 2024-T3 that bolts to a rib with two M6 bolts and carries an actuator's clevis on the sloping edge (not modelled). You will build it in FreeCAD, draw it, and then do what a checker does: verify, by an independent route, that the drawing is complete and consistent and that the model's mass properties are right. The independent route is the Python module of this lesson, which works only from the parameters and the dimension scheme, never from the CAD file.

Success criteria, written first: the CAD sketch reports "fully constrained" with no redundant constraints, and the Python constraint check agrees (0 DOF, 0 redundant); the CAD tool's volume, centre of mass and inertia agree with `mass_properties` to four significant figures; the drawing dimensions each feature once from the datum corner and carries units (mm), material, thickness and a general tolerance; the tolerance stack shows the bolt holes will assemble.

## 2. Understand and derive

### 2.1 What "fully constrained" means

A 2-D sketch is a set of points (and circle radii) with unknown coordinates. Each **point** has two degrees of freedom (DOF), each circle's radius one. The bracket has five outline corners and three hole centres, 16 coordinates, plus three radii: **19 unknowns**. Each dimension or geometric constraint (horizontal, vertical, coincident with the origin) is one equation. If the equations are independent and there are 19 of them, the sketch is **fully constrained**: one shape (locally). Fewer independent equations leave DOF (a corner can slide); more equations than independent ones are **redundant** (some dimension is implied by others), and if a redundant dimension's value disagrees with what the others imply, the sketch is **conflicting** and has no solution at all.

Formally, write the constraints as r(z) = 0 for the vector z of unknowns, and form the Jacobian J = ∂r/∂z (one row per constraint, one column per unknown). Then

**DOF = n_unknowns − rank J,  redundant = n_constraints − rank J.**   (1)

A CAD sketcher does exactly this: FreeCAD's Sketcher reports "Fully constrained", "Under constrained: k degrees of freedom" or "Redundant constraints" from the rank of its own Jacobian. The lesson's scheme is:

| # | constraint | equations |
|---|---|---|
| 1 | corner 0 at the origin (the datum) | 2 |
| 2 | bottom edge horizontal, overall width 60 | 2 |
| 3 | right edge vertical, height 15 | 2 |
| 4 | top corner of the slope at 20 across, 40 up, from the datum | 2 |
| 5 | left edge vertical, top edge horizontal | 2 |
| 6 | three hole centres, x and y from the datum | 6 |
| 7 | three hole radii (shown on the drawing as diameters ⌀6.6, ⌀6.6, ⌀12) | 3 |
| | total | 19 |

Every location is a **baseline dimension** from one datum (the lower-left corner, with the bottom and left edges as datum faces), which is what Ford's chapter 7 recommends for parts whose features must line up with a mating part. Adding the height of the sloping edge (25 mm, from y = 15 to y = 40) adds a **redundant** dimension: it is already 40 − 15. Removing the 20 mm leaves **one DOF**: the top corner can slide along the top edge and the slope changes. The code reports both.

### 2.2 Area properties by Green's theorem

For a region R bounded by a closed curve C traversed counter-clockwise, Green's theorem turns an area integral into a line integral: ∬_R (∂Q/∂x − ∂P/∂y) dA = ∮_C (P dx + Q dy). With Q = x, P = 0 the left side is the area A. On a straight edge from (x_i, y_i) to (x_{i+1}, y_{i+1}) the line integral is exact, and summing over the polygon's edges with c_i = x_i y_{i+1} − x_{i+1} y_i gives the **shoelace formulas**:

**A = ½ Σ c_i,  ∫x dA = (1/6) Σ (x_i + x_{i+1}) c_i,  ∫y dA = (1/6) Σ (y_i + y_{i+1}) c_i,**   (2)

**∫y² dA = (1/12) Σ (y_i² + y_i y_{i+1} + y_{i+1}²) c_i,  ∫xy dA = (1/24) Σ (x_i y_{i+1} + 2x_i y_i + 2x_{i+1} y_{i+1} + x_{i+1} y_i) c_i,**   (3)

and ∫x² dA like ∫y² dA in x (choose Q = x³/3 and so on). A clockwise polygon gives the same numbers with the opposite sign, so take the sign from A. The bracket's outline (0, 0), (60, 0), (60, 15), (20, 40), (0, 40) has A = 1900 mm². Holes are subtracted with the exact circle results: area πr², first moments πr²(x_c, y_c), second moments πr⁴/4 + πr² y_c² (parallel axis). The three holes remove 181.52 mm², leaving **A = 1718.48 mm²** with centroid **(25.618, 16.935) mm** and centroidal second moments I_xx = 207 561, I_yy = 472 268 and I_xy = −103 590 mm⁴. Principal values (`P05.05`, by `eigh`): 507 986 and 171 843 mm⁴, the major axis at 71.0° from x.

A CAD tool's B-rep model keeps circles exact. A mesh export (STL for a 3-D printer, or a finite-element mesh) replaces each hole by an n-gon. An inscribed n-gon has area (n/2)r² sin(2π/n) ≈ πr²(1 − 2π²/(3n²)), so each hole is too small by about 2π²/(3n²) of its area and the plate's net area comes out too large: by 1.05 % at n = 8 and 0.017 % at n = 64, falling as 1/n², which `test_polygon_circle_converges_to_exact_at_second_order` checks.

### 2.3 Mass properties of the extruded plate

Extrude the area by the thickness t along z. With density ρ: volume V = At, mass m = ρAt, centre of mass (x_c, y_c, t/2). The inertia tensor about the centre of mass (`P04.04`) of a prism of uniform section is

**J_xx = ρt I_xx,c + m t²/12,  J_yy = ρt I_yy,c + m t²/12,  J_zz = ρt (I_xx,c + I_yy,c),  J_xy = −ρt I_xy,c,**   (4)

and J_xz = J_yz = 0 by symmetry about the mid-plane. The t²/12 terms are the plate's own thickness; for a thin plate they vanish and J_zz = J_xx + J_yy (the perpendicular-axis theorem); for a thick one **J_zz = J_xx + J_yy − m t²/6**. For the bracket, t = 4 mm and ρ = 2780 kg/m³ (typical 2024-T3): **V = 6873.9 mm³, m = 19.11 g, centre of mass (25.618, 16.935, 2.000) mm**, J_xx = 2.3336, J_yy = 5.2771, J_zz = 7.5597 and J_xy = 1.1519 kg·mm².

The off-diagonal sign is a convention: (4) uses the **inertia tensor** convention J_xy = −∫xy dm. Some tools report the **product of inertia** +∫xy dm instead. Check which yours uses before comparing.

### 2.4 Tolerances and stack-ups

No dimension is made exactly. A **tolerance** ±t says the shop may make it anywhere in [nominal − t, nominal + t]. When a functional dimension y is a sum or difference of drawn dimensions, y = Σ a_i x_i, its variation is a **stack-up**:

**worst case: Δy = Σ |a_i| t_i,  root-sum-square: Δy = √Σ (a_i t_i)².**   (5)

Worst case guarantees every part assembles; RSS is the ±3σ band of y *if* each x_i is independent and normal with σ_i = t_i/3, and so lets 0.27 % of assemblies fall outside. If the real distributions are uniform over the tolerance band (a process that drifts across it), σ_i = t_i/√3 and the true 3σ band is √3 times the RSS figure. A seeded Monte Carlo makes the statistics concrete.

Two consequences for the drawing. First, **baseline versus chain**: the x distance between the lower bolt hole (10 mm from the datum) and the lightening hole (35 mm) is not dimensioned on the drawing; it is the difference of two baseline dimensions, each ±0.1 mm, so it varies by ±0.2 mm worst case (±0.141 RSS). If that distance mattered functionally it would be dimensioned directly, ±0.1. Dimension what the function needs. Second, the **floating-fastener rule**: when a bolt of diameter F passes through clearance holes of diameter H in two parts, each hole's axis may be off its true position by up to (H − F)/2 radially and the bolt still goes in. For M6 bolts (F = 6.0 mm) in 6.6 mm holes (the usual medium-fit clearance for M6), that is **0.30 mm**. The ±0.1 mm baseline tolerances on the hole coordinates put each hole within 0.1√2 = 0.14 mm of true position: the holes will assemble, with room.

### 2.5 Projection and the drawing

A multiview drawing shows the part by orthographic projection onto three perpendicular planes (Ford, chapters 1–2). In **third-angle projection** (standard in North America and in most aerospace drawing practice) the top view sits above the front view and the right-side view to the right of it; in **first-angle** (common in Europe and ISO practice) they are swapped. A small truncated-cone symbol in the title block says which. For a flat plate the front view (the plate face) carries almost every dimension; the edge view shows the thickness. Rules worth following (Ford, chapter 7): each feature dimensioned once, on the view that shows its true shape; dimensions from datums, not from other dimensions; holes by diameter with ⌀ and their centres located by centre lines; no dimension to hidden lines; units, material, thickness, general tolerance and projection symbol in the title block.

### 2.6 Where the models stop

- The mass properties are exact for a prism with through holes. Fillets, chamfers, countersinks and a non-uniform thickness break the formula; a CAD tool integrates the real solid and is the better reference then, provided the model is right.
- Full rank of the constraint Jacobian at the nominal shape is local: a set of distances can have a mirror-image solution far away (a corner that flips across an edge). CAD sketchers handle this by solving from your sketched starting point, which is why sketching roughly to shape matters.
- Linear stacks assume small variations and independent dimensions; tolerances on the same machine setup are correlated, and geometric tolerances (flatness, perpendicularity) add terms that coordinate tolerances do not capture.

## 3. Build the computational model

Open `starter/main.py`. The `Bracket` dataclass (outline, holes, thickness, density), `regular_polygon`, `point_in_polygon`, `bracket_sketch` (the lesson's dimension scheme as constraint tuples), `edge_distance_ratio` and `main` are given.

### Step 1: polygons and circles

`polygon_properties(P)` is (2)–(3) with `np.roll` for the next vertex; reverse the sign of everything if the area comes out negative. `circle_properties(xc, yc, r)` is the exact circle with the parallel-axis theorem.

### Step 2: the plate and its mass properties

`section_properties(b)` subtracts the holes (exact, or as n-gons for the convergence study) and returns the centroid and the centroidal moments; `mass_properties(b)` is (4). `raster_properties(b, h)` samples cell centres on an h-grid, keeps those inside the outline and outside every hole, and sums: a crude but completely independent route (no Green's theorem, no parallel axis). `principal_axes` uses `eigh` on the 2 × 2 area tensor.

### Step 3: tolerance stacks

`worst_case_stack` and `rss_stack` are (5); `monte_carlo_stack(nominal, sens, tol, n, seed, dist)` draws each dimension from a normal with σ = t/3 or a uniform on ±t with `numpy.random.default_rng(seed)`; `floating_fastener_position_tolerance` is (H − F)/2.

### Step 4: the dimension scheme

`constraint_residuals(z, cons, n_pts)` evaluates each constraint tuple (`fixx`, `fixy`, `dx`, `dy`, `dist`, `radius`, `horizontal`, `vertical`, `angle`); `constraint_jacobian` differentiates them by central differences; `constraint_status` applies (1) with the rank from the singular values (threshold 10⁻⁸ of the largest); `solve_sketch` rebuilds the geometry from a perturbed guess with `scipy.optimize.least_squares`, which is what a CAD sketcher's solver does.

Running `main()` prints

```
bracket outline area 1900.00 mm^2, holes 181.52 mm^2, net 1718.48 mm^2
centroid (25.618, 16.935) mm; raster h = 0.1 mm: (25.621, 16.934) mm, area 1718.24 mm^2
centroidal I_xx 207561, I_yy 472268, I_xy -103590 mm^4; raster 207529, 472199, -103576 mm^4
principal I1 507986, I2 171843 mm^4, I1 axis at 71.0 deg from x
volume 6873.9 mm^3, mass 19.11 g (2024-T3, 2780 kg/m^3), centre of mass (25.618, 16.935, 2.000) mm
inertia about the centre of mass [kg mm^2]: Jxx 2.3336, Jyy 5.2771, Jzz 7.5597, Jxy 1.1519
```

then the n-gon convergence, the edge-distance ratios (1.52, 1.52 for the bolt holes, 1.10 for the lightening hole), the tolerance and constraint checks of part 4, and writes `bracket_check.png` (the bracket drawn from its parameters with centre lines, three dimensions and the centroid: a check drawing to hold next to your CAD drawing).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-17-cad-technical-drawings`. The tests:

1. **Rectangle**: area, centroid, bh³/3 and hb³/3 about the edges, b²h²/4.
2. **Triangle**: area, centroid at a third, bh³/12 about the base; independence of orientation and starting vertex.
3. **Circle**: n-gon area error falls at order 2.00 ± 0.05; centroid exact at every n.
4. **Bracket against raster** (independent check): area to 10⁻³, centroid to 0.01 mm, second moments to 2 × 10⁻³; a 256-gon hole agrees with the exact one to 10⁻⁴.
5. **Holes**: net area = outline − Σπr²; a hole near the corner moves the centroid away from it.
6. **Principal axes**: I₁ + I₂ = I_xx + I_yy; the moment about the computed axis equals I₁.
7. **Mass properties**: V = At, m = ρV, J symmetric positive definite, J_zz = J_xx + J_yy − mt²/6, doubling t doubles m and moves the centre of mass to t/2.
8. **Box formula** for a solid rectangular plate, and zero product of inertia.
9. **Stacks**: worst case ≥ RSS, equal for one dimension; Monte Carlo 3σ = RSS to 1 %, coverage 0.9973 ± 0.001, uniform bounded by the worst case with σ = RSS/√3.
10. **Floating fastener** (H − F)/2.
11. **Scheme**: residuals zero at nominal; 0 DOF, 0 redundant; one extra dimension is redundant; dropping one leaves 1 DOF; dropping the x datum leaves the part free to slide.
12. **Rebuild**: `solve_sketch` from a perturbed guess returns the nominal geometry to 10⁻⁸ mm, and the rebuilt outline has the right area.
13. **Distance and angle** constraints on a 3-4-5 triangle.
14. **main** writes the figure.

### Compare

The established tool here is the CAD program itself, and part 4 is the comparison. FreeCAD's solid modeller integrates the exact B-rep solid, so its volume and centre of mass should match `mass_properties` to at least five significant figures for this prism; a difference in the third or fourth figure means the model and the parameters disagree (a hole of the wrong size, a pad of the wrong length), which is exactly the kind of error the check exists to catch.

## 4. Experiment at home, or investigate real data

**Mode: design-build** in free CAD software, with a computational check. FreeCAD menu names below are as in FreeCAD 1.0; earlier versions put some commands one menu level higher, and some names differ slightly between versions ("Create polyline" or "Polyline"), so look for the same command and icon. In another tool, do the same steps with its equivalents.

**Predict first.** Before opening the CAD tool, write down from part 2: the volume (6873.9 mm³), the mass (19.11 g) and the centre of mass (25.62, 16.94, 2.00) mm you expect it to report, and how many constraints the sketch will need (19, counting the origin lock as two).

**Model.**

1. Start FreeCAD, create a new document (**File → New**), and switch to the **Part Design** workbench (workbench selector in the toolbar). Create a body (**Part Design → Create body**).
2. Optional but recommended: create a spreadsheet (**Spreadsheet** workbench, **Spreadsheet → Create spreadsheet**) with named cells `width = 60`, `height = 40`, `right = 15`, `slope_x = 20`, `t = 4`, `bolt_d = 6.6`, `light_d = 12` and the hole coordinates, so that every dimension below is an expression such as `=Spreadsheet.width` (FreeCAD manual, "Using spreadsheets").
3. **Part Design → Create sketch**, on the XY plane. Draw the outline roughly to shape with **Sketch → Sketcher geometries → Create polyline**, closing it on the first point, and three circles with **Create circle**.
4. Constrain: put the first corner on the origin (**Constrain coincident** with the sketch's origin point); bottom and top edges **Constrain horizontal**, left and right edges **Constrain vertical**; then the dimensions with **Constrain horizontal distance**, **Constrain vertical distance** and **Constrain diameter**, exactly as the table in 2.1, every location measured from the origin corner. Watch the Sketcher's solver messages: it must say **Fully constrained** and must not list redundant or conflicting constraints. Close the sketch.
5. **Part Design → Create an additive feature → Pad** the sketch by 4 mm (the holes, being closed circles inside the outline, become through holes automatically; alternatively draw only the outline, pad it, and add the holes with **Pocket** or **Hole**).
6. Read the mass properties. Open the Python console (**View → Panels → Python console**) and type, with your body's name:
   ```
   s = App.ActiveDocument.getObject("Body").Shape
   s.Volume, s.CenterOfMass, s.MatrixOfInertia
   ```
   The volume is in mm³, the centre of mass in mm, and the matrix of inertia is for unit density (mm⁵): multiply by ρ = 2.78 × 10⁻⁶ kg/mm³ to compare with J in kg·mm². Record all three.

**Draw.**

7. Switch to the **TechDraw** workbench. **TechDraw → Pages → Insert Default Page** (A4 landscape is enough). Select the body and insert a projection group (**TechDraw → Views → Insert Projection Group**) with the front (plate face), top and right views; set the page's projection type to third angle in its properties.
8. Dimension the front view from the datum corner with **TechDraw → Dimensions → Insert Horizontal Dimension**, **Insert Vertical Dimension** and **Insert Diameter Dimension**; give the thickness (4) on an edge view. Add centre lines to the holes. Do not dimension anything twice: compare your drawing's dimensions one by one with the table in 2.1.
9. Fill the title block (double-click its editable fields): part name, your name, date, units mm, material 2024-T3 sheet 4.0 thick, general tolerance ±0.1 mm, third-angle projection. Export the page as PDF (**File → Export** with the page selected, or print to PDF).

**Check.**

10. Enter your drawing's dimension scheme as constraint tuples (start from `bracket_sketch` and change it if your scheme differs) and run `constraint_status`. It must report 0 DOF and 0 redundant; if the CAD sketcher and the Python check disagree, find out why.
11. Compare the CAD tool's volume, centre of mass and inertia with `mass_properties`.
12. Run the tolerance stack of `main` for the assembly gap (a 62.0 ± 0.10 mm housing slot, the 60.0 ± 0.05 mm bracket, two 0.50 ± 0.02 mm shims).

Raw data format: `cad_check.csv` with columns `quantity,cad_value,python_value,unit,relative_difference` (volume, x_c, y_c, z_c, J_xx, J_yy, J_zz, J_xy); `scheme.txt` with your constraint tuples and the `constraint_status` output; the drawing PDF.

What you should find: the CAD and Python values agree to five or more significant figures (V = 6873.9 mm³, centre of mass (25.618, 16.935, 2.000) mm), and J_xy appears with the sign convention of your tool. The constraint check reports 19 unknowns, 19 constraints, rank 19. The stack gives a nominal gap of 1.000 mm with ±0.190 mm worst case and ±0.115 mm RSS; a seeded Monte Carlo (seed 17, 200 000 samples) gives 3σ = 0.1154 mm and 99.73 % of gaps inside the RSS band; with uniform distributions only 92.6 % are, and the 3σ figure (0.2003 mm) even exceeds the worst-case band, which shows that 3σ is not a bound for a uniform distribution. The hole-to-hole x spacing varies ±0.200 mm worst case under baseline dimensioning, against ±0.100 mm if dimensioned directly.

## 5. Validate and improve

| quantity | theory (§2) | code / CAD | difference | cause |
|---|---|---|---|---|
| net area | 1900 − 181.52 = 1718.48 mm² | raster 1718.24 mm² | 1.4 × 10⁻⁴ | boundary cells, h = 0.1 mm |
| centroid | (25.618, 16.935) mm | raster (25.621, 16.934) mm | 0.003 mm | same |
| I_xx,c, I_yy,c, I_xy,c | 207 561, 472 268, −103 590 mm⁴ | raster 207 529, 472 199, −103 576 | ≤ 1.6 × 10⁻⁴ | same |
| hole as 64-gon | exact circle | net area +1.70 × 10⁻⁴ | 1/n² | inscribed polygon |
| J_zz | J_xx + J_yy − mt²/6 | 7.5597 kg·mm² | < 10⁻¹² | — |
| Monte Carlo 3σ (normal) | RSS = 0.115 mm | 0.1154 mm | 0.3 % | sampling, 2 × 10⁵ |
| scheme | 19 independent equations | rank 19 | — | — |
| CAD volume, centre of mass | 6873.9 mm³, (25.618, 16.935, 2.000) | your FreeCAD values | expect < 10⁻⁵ | — |

The raster integration is the independent check of the formulas, the n-gon study of the holes, the Monte Carlo of the stack statistics, and your CAD model of all of them together. A mass from CAD that disagrees with a hand formula is not "CAD error": it is a modelling mismatch to be found, and the most common causes are a pad of the wrong length, a hole drawn as a radius where a diameter was meant, and a density left at the tool's default.

What to improve: drive the CAD model from a spreadsheet (Extension 1) so that drawing, model and check share one set of parameters; dimension the holes with a position tolerance (Extension 3). What not to improve: the raster resolution, which only confirms what the exact formulas already give.

## 6. Defend your work

Submit: the CAD file and the drawing PDF; `cad_check.csv`; `scheme.txt`; the tolerance stack table; the completed `main.py` with passing tests and `bracket_check.png`; and a half-page drawing-check note (template in `projects/_template/report-template.md`, sections 1, 4, 5, 6 and 7) stating whether the drawing is complete and consistent, the mass properties and their agreement, the stack-up result and the method used, and one thing you would change on the drawing.

### Adversarial self-checks

- **What if the thickness doubles?** Mass and volume double, the centre of mass moves to z = 4 mm, J_zz doubles exactly, J_xx and J_yy more than double (the t² term grows fourfold). The constraint scheme does not change: thickness is a pad length, not a sketch dimension, and it must still appear once on the drawing.
- **What if a tolerance goes to zero?** It drops out of both stacks; but a zero tolerance cannot be made or inspected, so a drawing with an untoleranced dimension relies on the general tolerance in the title block. Check that there is one.
- **Which assumption fails first in the real world?** The plain prism: a real bracket would have corner fillets, possibly a countersink and edge breaks. The CAD tool then becomes the reference and the formulas the sanity check (fillets change the mass by well under 1 % here).
- **What measurement could falsify the model?** Weigh the made part: 19.11 g predicted; a kitchen scale with 0.1 g resolution resolves 0.5 %, enough to catch a wrong thickness (3.8 mm stock would give 18.15 g) or a missing hole.
- **Does the result survive a change of method?** Green's theorem, raster sampling and the CAD solid modeller are three different integrations; the constraint status from the Jacobian rank and from the CAD sketcher are two different solvers.
- **Which decision changes if the shop's process is uniform rather than normal?** The RSS band covers only 92.6 % of assemblies instead of 99.73 %; if the gap must never close, design to worst case or tighten the shim tolerance.

## Check yourself

1. How many degrees of freedom does a sketch of five points and three circles have before any constraint is added, and how many independent constraints make it fully constrained?
2. What is the difference between a redundant and a conflicting dimension?
3. Write the shoelace formula for the area of a polygon. What happens to its sign if the vertices are listed clockwise?
4. A plate's area is 1718.48 mm², its thickness 4 mm and its density 2780 kg/m³. What is its mass?
5. Why is J_zz of a thick plate smaller than J_xx + J_yy, and by how much?
6. Three dimensions with tolerances ±0.10, ±0.05 and ±0.02 mm add up to a gap. Give the worst-case and RSS variations.
7. What positional error can each hole have when 6.0 mm bolts pass through 6.6 mm clearance holes in both parts?
8. In third-angle projection, where does the top view sit relative to the front view?
