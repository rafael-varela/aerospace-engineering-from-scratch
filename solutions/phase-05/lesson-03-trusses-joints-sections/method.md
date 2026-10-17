# P05.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Apex a = (1 + s/30, 3). 2. ΣM_pin: 4 R₁ = a_x × 1000; vertical sum for the pin reaction. 3. Joint 1 (roller): vertical equilibrium gives the force in member 1–2 (direction (a − joint 1)/|…|); horizontal gives member 0–1. 4. Joint 2 or joint 0 for member 0–2. 5. `Truss(...)`, `solve_joints`, `joint_residual`.

## 2.
(a) Two non-collinear members at an unloaded joint: resolve along the normal to one member, only the other has a component. Three members, two collinear: resolve normal to their line. (b) Move a diagonal: same m, one panel double-braced (indeterminate), another empty (mechanism). Rank deficit counts both. (c) Free body left of the cut: moments about the top end of the cut diagonal leave only the bottom chord (lever h); about its bottom end, only the top chord; vertical sum, only the diagonal.

## 3.
1. m + r = 18 = 2j. 2. Build nodes, members, fixed (pin at 0, roller at 4), loads at upper joints. 3. `solve_joints`; argmin/argmax. 4. Section: cut = indices of the three members, side = joints left of the cut (lower 0, 1, 2 and upper 5, 6); `external_reactions_2d`; `solve_section`. 5. V at the cut = reaction − loads to the left. 6. M_max at mid-span; divide by h.

## 4.
1. span = 4 + s/10. 2. Sweep; argmin of the volume index; the matching chord force. 3. At 0.25 span: group extremes. 4. A = `tube_area(0.01905, 0.000889)`, I = `tube_second_moment(...)`; σ = T/A; `euler_load_pinned(200e9, I, span/6)`; compare.

## 5.
(a) `fe, Re = solve_joints(em)`; extremes; `member_geometry(em)[0]` for lengths; stresses; Euler load. (b) New `Truss` with the extra member; `determinacy`; `stiffness_solve(em2, E*A)` and with the last EA × 4. (c) `np.delete(em.members, k, axis=0)`; `determinacy`.

## 6.
(a) Slender members, joints at centroid intersections; misses secondary bending. (b) The bearer is a beam-column; its end forces have transverse components; the truss joints at its ends see extra loads. (c) Buckling of the longest compression member. (d) Indeterminate; load shares by stiffness; thermal and fit-up stresses; a stiff new path can overload its neighbours or reverse forces.
