# P04.09  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The force is qvB because v ⊥ B; equate it to mv²/r for the radius; the period is 2πm/(qB) and does not depend on v. The right-hand rule gives the direction of v × B; the electron's charge flips it.
2. (a) On the axis every element is the same distance from the field point, and for each element there is an opposite one whose transverse contribution cancels. (b) A slice dz′ of the solenoid is a loop carrying nI dz′; substitute z − z′ = R tan u. (c) The sides of length a are always perpendicular to B and carry equal and opposite forces IaB that are not on the same line; the sides of length b carry forces along the rotation axis that cancel on one line.
3. The motional emf is vBL; the current flows through the 100 Ω loop; the force on the current-carrying tether is ILB. Lenz: the induced current's force must oppose the motion that induced it.
4. `field_along_orbit(400e3, 90.0)` and `np.linalg.norm(..., axis=1)`. In the table set I = m_req/(N side²). The saturation rate m_max/(kB) tells you whether the coil saturates from 5°/s.
5. `loop_points(R, 720, center=[0, 0, d/2])` for one loop and `center=[0, 0, -d/2]` for the other; add the two fields. For (b), expand B(z) of each loop about the midpoint; the odd terms cancel by symmetry.
6. Put numbers on each claim: (R_E/r)³ for GEO; μ₀m/(2πz³) at z = 0.1 m; τ·B = 0 for any combination of coils.
