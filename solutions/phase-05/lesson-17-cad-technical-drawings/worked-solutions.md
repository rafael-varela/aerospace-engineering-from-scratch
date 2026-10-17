# P05.17  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. An angle section (s = 15: H = 33 mm)

**Assumptions.** Sharp corners (no root fillet), units mm.

(a) Rectangles 40 × 4 (A₁ = 160, centroid (20, 2)) and 4 × 29 (A₂ = 116, centroid (2, 18.5)). **A = 276 mm²; x̄ = (3200 + 232)/276 = 12.435 mm; ȳ = (320 + 2146)/276 = 8.935 mm.**

(b) `polygon_properties` gives the same area and centroid, and **I_xx,c = 26 651 mm⁴, I_yy,c = 43 276 mm⁴, I_xy,c = −19 972 mm⁴.**

(c) **I₁ = 56 596 mm⁴, I₂ = 13 330 mm⁴**, the major axis at **56.3°** from x. I_xy ≠ 0 because the section has no axis of symmetry through its centroid parallel to x or y: the material in the horizontal leg lies at x > x̄ with y < ȳ and that in the vertical leg at x < x̄ with y > ȳ, so ∫(x − x̄)(y − ȳ) dA is negative.

**Check.** I₁ + I₂ = 69 926 = I_xx,c + I_yy,c. The composite-rectangle sums (each rectangle's own I plus its parallel-axis term) give the same three centroidal values. An equal-leg angle would have its principal axes at exactly 45°; here the eigenvector gives 56.3°, and the trustworthy check is the invariant sum and the composite calculation, not an intuition about which leg the axis should lean towards.

**Common mistakes.** Using the vertical leg's full height H for A₂ (double-counting the 4 × 4 corner); forgetting the parallel-axis terms when moving I to the centroid.

## 2. Derivations

(a) With Q = x, P = 0, Green gives A = ∮x dy. On edge i, x = x_i + τΔx_i, dy = Δy_i dτ: ∫₀¹ (x_i + τΔx_i)Δy_i dτ = ½(x_i + x_{i+1})(y_{i+1} − y_i). Summing, the terms ½(x_i y_{i+1} − x_{i+1} y_i) survive and the others (x_i y_i, x_{i+1} y_{i+1}) cancel in pairs around the closed polygon: **A = ½Σ(x_i y_{i+1} − x_{i+1} y_i)**. With Q = x²/2: ∫x dA = ½∮x² dy = ½Σ∫₀¹(x_i + τΔx_i)²Δy_i dτ = (1/6)Σ(x_i² + x_ix_{i+1} + x_{i+1}²)Δy_i, which rearranges (using the closed sum again) to **(1/6)Σ(x_i + x_{i+1})(x_i y_{i+1} − x_{i+1} y_i)**.

(b) About the centre of mass, with z from −t/2 to t/2: J_xx = ∫(y² + z²) dm = ρt I_xx,c + ρA∫z² dz = ρt I_xx,c + m t²/12; J_yy = ρt I_yy,c + m t²/12; J_zz = ∫(x² + y²) dm = ρt(I_xx,c + I_yy,c). Hence **J_zz = J_xx + J_yy − m t²/6**. The perpendicular-axis theorem J_zz = J_xx + J_yy is the limit t → 0, where all the mass lies in the plane z = 0.

(c) For independent x_i, Var(Σa_ix_i) = Σa_i²σ_i². Normal with σ_i = t_i/3: σ = √Σ(a_it_i)²/3, so **3σ = √Σ(a_it_i)² = RSS**. Uniform on ±t_i: σ_i = t_i/√3, so 3σ = 3√Σa_i²t_i²/3 = **√3 × RSS**.

**Common mistakes.** Signs from the orientation (a clockwise polygon gives −A); applying the parallel-axis theorem in the wrong direction.

## 3. A slotted bracket

**Assumptions.** Typical 7075-T6 density 2810 kg/m³; t = 3 mm; slot fully inside the plate.

(a) **A_slot = 14 × 6.6 + π × 3.3² = 92.40 + 34.21 = 126.612 mm²**; the 400-vertex polygon gives 126.611 mm² (inscribed arcs, 10⁻⁵ low).

(b) Plate 2400 mm², net **A = 2273.39 mm²**; the slot's centroid is at its centre (50, 15) by symmetry, so x̄ = (2400 × 40 − 126.61 × 50)/2273.39 = **39.443 mm**, ȳ = **15.000 mm**. Mass **ρAt = 19.16 g**, against 20.23 g for the solid plate: the slot saves **1.07 g** (5.3 %).

(c) The slot is two arcs and two lines: 4 + 4 × 2 endpoints, but with the arcs' centres and radii as unknowns; FreeCAD's slot tool applies tangency and equal-radius constraints automatically. To fix it you need: one radius (or the width), the centre-to-centre distance (or overall length), the slot's orientation (horizontal) and the position of its centre (two dimensions from the datum): **five independent dimensions and constraints** beyond the automatic ones. On the drawing: "2 × R3.3" or "6.6 wide", "14" between arc centres (or overall length 20.6), and the centre located from the datum (50, 15); the overall length then must not also be given.

**Common mistakes.** Counting two semicircles as two circles; locating the slot from its end instead of its centre and then also giving the centre.

## 4. Changing the dimension scheme

(a) The signed angle at (60, 15) from the line to (60, 0) to the line to (20, 40) is **−2.1294 rad (122.0° in magnitude)**, so the constraint is ("angle", 1, 2, 3, −2.1294). With it replacing ("dx", 0, 3, 20.0): **19 unknowns, 19 constraints, rank 19, 0 DOF, 0 redundant.** The angle is as good a dimension as the 20 mm; it suits a part whose sloping edge is cut at a set angle, and it puts the 20 mm out of direct control.

(b) With ("dy", 5, 6, 20.0) instead of ("dy", 0, 6, 30.0): still **0 DOF, 0 redundant** (a chain dimension still fixes the hole). The upper hole's height above the datum is now 10 + 20, with **±0.2 mm** worst case (±0.141 RSS) instead of ±0.1; its distance from the lower hole is now held to ±0.1 instead of ±0.2. Which is right depends on which distance the mating part cares about: for two bolts into a rib's hole pair, the hole-to-hole distance is the functional one, and the chain is the better choice.

(c) **20 constraints, rank 19: 0 DOF, 1 redundant.** The value 26 contradicts the 25 implied by 40 and 15: the system is conflicting. `solve_sketch` finds the least-squares compromise and leaves a **maximum residual of 0.25 mm**, spread over the conflicting dimensions. FreeCAD's Sketcher reports "conflicting constraints" and names them, and refuses to solve until one is removed. A non-zero residual after a solve is how the Python check tells you the same thing.

**Common mistakes.** Taking the angle with the wrong sign or from the wrong pair of lines, which makes the scheme inconsistent; deleting a constraint at random to fix a conflict instead of finding which dimension is wrong.

## 5. Tolerance stack with your own shim (s = 15: shim ±0.095 mm)

**Assumptions.** Linear stack, independent dimensions.

(a) Gap = 49.5 − 48.0 − 1.0 = **0.500 mm**. Worst case ±(0.10 + 0.05 + 0.095) = **±0.245 mm**, smallest gap **0.255 mm** (fails 0.30). RSS ±√(0.01 + 0.0025 + 0.009025) = **±0.147 mm**, smallest gap **0.353 mm** (passes).

(b) Normal, seed 15: 3σ = **0.1471 mm**, fraction below 0.30 mm **1 × 10⁻⁵** (2 of 200 000). Uniform: 3σ = **0.2541 mm**, fraction below 0.30 mm **0.40 %**.

(c) Worst case needs 0.10 + 0.05 + t ≤ 0.20, so **t ≤ 0.05 mm** on the shim. For accepting the original design: the RSS figure holds if the shop's processes are centred and roughly normal, and then failures are one in a hundred thousand. Against: the requirement says "never"; the uniform case, which is what a drifting process looks like, gives 0.4 % of assemblies too tight; a shim is cheap to make to ±0.05. Tighten the shim tolerance; it costs almost nothing and removes the dependence on the distribution.

**Common mistakes.** Adding the tolerances of quantities that subtract with a minus sign (tolerances always add in the worst case); quoting an RSS result without the distribution assumption.

## 6. Judgment: a drawing review

A model review, in outline. **Findings.** (1) The sloping-edge height 25 is redundant with 40 and 15: an inspector measuring a part with height 40.08 and right edge 14.95 will find the 25 out of tolerance even if both others are in; delete it (or mark it as reference, "(25)"). (2) The lightening hole is chained from the lower bolt hole (25, 5); this is complete but accumulates the bolt hole's tolerance; locate it from the datum at (35, 15) unless its position relative to the bolt hole is functional (it is not; it only removes mass). (3) The top-edge length 20 is acceptable: the left edge is a datum face, so it is the same as 20 from the datum. (4) No general tolerance: every unmarked dimension is uninspectable; add "unless otherwise stated ±0.1 mm". (5) No material, thickness note or units in the title block: add "2024-T3 sheet, 4.0 thick, dimensions in mm". (6) First-angle symbol with the top view above the front view: the symbol and the layout contradict each other; in first angle the top view goes below. A shop could make a mirror-image part from a misread projection. Use third angle and its symbol. **Corrected scheme** = the lesson's table of 2.1. **Check.** The colleague's scheme gives 19 unknowns, 20 constraints, rank 19: **1 redundant**; with the 25 removed, 0 DOF and 0 redundant (the chained hole location does not change completeness, only tolerance accumulation).

## Check-yourself answers

1. 5 × 2 + 3 × 3 = **19 DOF**; **19 independent constraints** make it fully constrained (the origin lock counts as two).
2. A redundant dimension is implied by others (the rank does not increase); if its value agrees with what the others imply it is merely redundant, if it disagrees the sketch is **conflicting** and has no solution.
3. **A = ½Σ(x_i y_{i+1} − x_{i+1} y_i)**, summed cyclically. Clockwise order gives **−A**.
4. m = 2780 kg/m³ × 1718.48 × 10⁻⁶ m² × 0.004 m = **19.11 g**.
5. J_xx and J_yy each include the thickness term m t²/12 about the mid-plane, J_zz does not: **J_zz = J_xx + J_yy − m t²/6**.
6. Worst case **±0.17 mm**; RSS √(0.01 + 0.0025 + 0.0004) = **±0.114 mm**.
7. **(6.6 − 6.0)/2 = 0.30 mm** radially (a ⌀0.60 mm position zone) for each part.
8. **Above** it (and the right-side view to the right of the front view).
