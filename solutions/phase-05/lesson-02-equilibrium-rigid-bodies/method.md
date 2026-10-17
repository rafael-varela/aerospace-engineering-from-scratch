# P05.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. P = 1000 + 100s. 2. Unknown counts: (a) 3, (b) 3, (c) 4, (d) 3, (e) 3, (f) 3. 3. Geometry: (d) all parallel; (f) the link's line passes through the pin (concurrent). 4. Build the supports with `pin`, `roller`, `fixed_2d`, `link`; `classify_supports(equilibrium_matrix(U, dim=2))`. 5. Determinate cases: moments about one support point, then vertical and horizontal sums; check with `solve_reactions(U, [[1.5, 0]], [[0, -P]], dim=2)`.

## 2.
(a) Every column has moment (r − O) × u = 0 about the meeting point O, so the moment row of A is a combination of the others: rank ≤ 2. Free motion: rotation about O. Self-stress: the three link forces in proportions that sum to zero (a closed force triangle). (b) dim null(A) = n − ρ (self-stress states); dim null(Aᵀ) = n_eq − ρ (loads orthogonal to the range of A cannot be balanced). (c) ΣM_root: T × (y_s h/ℓ_s) = L ȳ − nW s/2. Horizontal pin force = T y_s/ℓ_s = (L ȳ − nW s/2)/h.

## 3.
1. 3 + 2 + 1 = 6 unknowns; check rank. 2. Moment about the x axis of each load: M_x = y F_z − z F_y. 3. Strut direction u = (frame − door)/|frame − door|; unit moment y u_z − z u_y at the door end. 4. T = −ΣM_x(loads)/(unit moment); sign gives the sense. 5. `solve_reactions` with `ball_socket`, `hinge(..., axial=False)`, `link`. 6. `residual`.

## 4.
1. m = 900 + 20s; L = ½ × 3.8 × m g. 2. `strut_wing_reactions` and `strut_force_closed_form` at y_s = 2.6. 3. Sweep `np.linspace(1.0, 4.6, 361)`; resultant `np.hypot(root_x, root_y)`; argmin. 4. Negative case: lift_half = −0.4 L, n = −1.52. 5. Reasons: length, buckling, wing torsion and aileron/flap clearance, drag.

## 5.
(a) `solve_reactions(U, np.zeros((1, 3)), np.zeros((1, 3)), np.array([[450.0, 0, 0]]))`; max |force|; × sin(spread); `classify_supports(...)["cond"]`. (b) `A = equilibrium_matrix(U7)`, `b = load_vector(*loads)`, `x = lstsq(A, -b)`; check A x + b; null vector from the SVD. The real share depends on the links' stiffnesses (P05.03 stiffness method). (c) `hexapod_mount(drop_link=k)` for k = 0..5.

## 6.
(a) Determinacy: no stress from thermal or structural distortion reaches the instrument. (b) A second load path after a link failure; requires a stiffness (compatibility) analysis and a check of the remaining structure with one link gone. (c) It cannot resist a horizontal load at all, however small; and it is also indeterminate for vertical loads. (d) Forces of order condition number × load/lever; small geometric errors change the forces in proportion.
