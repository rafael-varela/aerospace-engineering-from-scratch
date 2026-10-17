# P05.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an angle section by hand and by code)

An extruded angle (L-section) stringer has the outline (0, 0), (40, 0), (40, 4), (4, 4), (4, H), (0, H) in mm, with H = 30 + s/5.

(a) By splitting it into two rectangles, find the area and the centroid.
(b) With `polygon_properties`, find the area, the centroid and the centroidal I_xx, I_yy, I_xy. Confirm (a).
(c) Find the principal second moments and the angle of the major principal axis with `principal_axes`. Why is I_xy not zero for this section?

Answer form: area and centroid twice, three second moments, two principal values and an angle.

## 2. (Derivation: Green's theorem and the thick plate)

(a) From Green's theorem with suitable P and Q, derive the shoelace formulas for the area and for ∫x dA of a polygon, by evaluating the line integral exactly on one straight edge.
(b) Show that for a prism of thickness t and uniform density, J_zz = J_xx + J_yy − m t²/6 about the centre of mass, and say why the perpendicular-axis theorem J_zz = J_xx + J_yy holds only for a thin plate.
(c) Show that if each of n independent dimensions is normal with σ_i = t_i/3, the sum Σ a_i x_i has a ±3σ band of √Σ(a_i t_i)², and that if instead each is uniform on ±t_i, the 3σ figure is √3 times larger.

Answer form: three derivations.

## 3. (Unfamiliar scenario: a slotted bracket)

An access-panel bracket is an 80 × 30 mm plate of 7075-T6 (typical density 2810 kg/m³), 3 mm thick, with one obround slot: width 6.6 mm, centre-to-centre length 14 mm, centred at (50, 15) mm from the lower-left corner. A slot lets the bracket be adjusted along x during assembly.

(a) Compute the slot's area exactly (a rectangle plus two semicircles) and check it against a polygon approximation of the slot with 400 vertices in `polygon_properties`.
(b) Compute the net area, the centroid and the mass of the bracket. How much mass does the slot save?
(c) In FreeCAD, which sketch constraints define the slot (count them), and how would you dimension it on the drawing so that it is located once and sized once?

Answer form: two slot areas, net area, centroid, two masses, a constraint count and a dimensioning scheme.

## 4. (Code: changing the dimension scheme)

Start from `bracket_sketch(Bracket())`.

(a) Replace the 20 mm horizontal location of the slope's top corner by an angle constraint at the corner (60, 15) between the right edge and the sloping edge. Compute the angle the constraint must take, and show with `constraint_status` that the new scheme is still fully constrained without redundancy.
(b) Replace the upper bolt hole's vertical baseline dimension (30 mm from the datum) by a chain dimension of 20 mm from the lower bolt hole. Is the scheme still complete? What changes in the tolerance of the upper hole's height above the datum if both dimensions are ±0.1 mm?
(c) Add the height of the sloping edge as 26 mm (the drawing's other dimensions imply 25). What does `constraint_status` report, and what does `solve_sketch` do? How does a CAD sketcher report the same situation?

Answer form: an angle, three constraint-status reports, two tolerances, and the residual of the least-squares solve.

## 5. (Tolerance stack with your own shim)

A bracket of width 48.0 ± 0.05 mm sits in a housing slot of 49.5 ± 0.10 mm with one shim of 1.00 ± (0.02 + 0.005s) mm. The remaining gap must never be less than 0.30 mm.

(a) Compute the nominal gap, the worst-case and RSS variations, and the smallest gap by each method.
(b) Run `monte_carlo_stack` (200 000 samples, seed s) with normal and with uniform distributions; report the 3σ value and the fraction of assemblies with a gap below 0.30 mm.
(c) Which shim tolerance would make the worst case meet the requirement? Is the RSS result enough to accept the original design? Argue both sides.

Answer form: three numbers, two pairs of numbers, one tolerance and a short argument.

## 6. (Judgment: a drawing review)

A colleague's drawing of the lesson's bracket carries these dimensions: overall width 60; overall height 40; right-edge height 15; sloping-edge height 25; top-edge length 20; bolt holes ⌀6.6 at (10, 10) and (10, 30) from the lower-left corner; lightening hole ⌀12 located 25 mm to the right of the lower bolt hole and 5 mm above it; thickness 4; no general tolerance; first-angle projection symbol, but the top view drawn above the front view.

Write a review (about 250 words): list every problem (missing, redundant, chained or ambiguous dimensions, missing notes, projection inconsistency), say which ones could cause a wrong or rejected part and how, and give the corrected dimension list. Check your corrected list with `constraint_status`.

Answer form: a list of findings, a corrected dimension list and its constraint-status output.
