# P05.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. h = 40 + s. 2. Flange: A = 60 × 8, ȳ = h − 4, I_c = 60 × 8³/12. Web: A = 6(h − 8), ȳ = (h − 8)/2, I_c = 6(h − 8)³/12. 3. ȳ = ΣAȳ/ΣA. 4. Ix = Σ(I_c + A d²). 5. Iy = 8 × 60³/12 + (h − 8) × 6³/12. 6. S = Ix/c for c to the top and bottom; r_x = √(Ix/A). 7. `composite_section([rect(...), rect(...)])`; polygon of eight vertices; `polygon_section`.

## 2.
(a) On the edge, dx = Δx dτ, dy = Δy dτ; x dy − y dx = (x_i Δy − y_i Δx) dτ = c_i dτ; integrate. For Q_y: the triangle (origin, i, i + 1) has signed area c_i/2 and centroid at one third of the sum of its vertices; Q_y = Σ (area × centroid x). (b) y = y_c + d; ∫ y_c dA = 0 by the definition of the centroid. For Ixy, x = x_c + d_x, y = y_c + d_y. (c) u = x cos θ + y sin θ, v = −x sin θ + y cos θ; expand ∫v² dA, ∫u² dA, ∫uv dA; add and multiply for the invariants; set I_uv = 0. Normal n = (n_x, n_y): distance from the axis is n·r; square and integrate.

## 3.
1. Outer arc: angles π/2 to 3π/2 at R = 30 mm (the polygon closes along x = 0). 2. Inner arc at 27.5 mm between ±arccos(−2/27.5) going round through π. 3. `polygon_section(outer, [inner])`; repeat with 50 and 200 points. 4. Thin-walled: r_m = 28.75 mm; arc area π r_m t, centroid −2r_m/π; web area 2r_m × 2, centroid −1 mm. 5. Masses ρA; combined x = Σ m x/Σ m; add 30 mm for the distance from the nose. 6. A dense mass far forward moves the centre of mass more per kilogram than material spread round the nose, which also sits partly aft.

## 4.
1. h = 20 + s/3 mm. 2. Three routes; percentage errors of the thin-walled values. 3. `principal_axes`, `principal_axes_eig`; the eigenvector of the larger eigenvalue is the normal of the I_max axis, so the axis is at the eigenvector's angle minus 90°. 4. `rotate_axes(..., theta_p)`. 5. r_min = √(I_min/A). 6. `minimize_scalar(lambda b: -ratio(b), bounds=(0.004, 0.04), method="bounded")`.

## 5.
1. t = 1.0 + s/30 mm. 2. Three idealisations; errors relative to full; ρ × A. 3. Skins: two `rect(0.4 − 0.002, t, 0, ±(0.12 − t)/2)`; their I_c + A d²; divide by full Ix. 4. Loop over n_stringers.

## 6.
1. Polygon Ix at t ± 0.1 mm against the thin-walled error. 2. `i_section_parts` plus two `rect(0.0048, 0.006, ±0.015, (0.15 − 0.006)/2, sign=-1)`; centroid shift and Ix change. 3. Fillets add material at the web–flange junction, close to the flange; small. 4. Preliminary: nominal or minimum thickness; stress check: net section at the holes.
