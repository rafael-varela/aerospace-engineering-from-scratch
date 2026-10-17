# P04.04  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `inertia_solid_cylinder(m, 0.1, 0.02)`; h = I₃ω, T = ½I₃ω²; `parallel_axis(I, m, [0.5, 0, 0])`; for (d) the total angular momentum about axis 3 starts at zero.
2. (a) Apply (4) to L = Iω with I constant in body axes. (b) Dot Euler's equations with ω for the energy, with Iω for |L|². (c) Keep only first-order terms in ε₁, ε₃; differentiate one equation and substitute the other. (d) Steady motion: dL/dt is the precession Ω ẑ × L. (e) At a turning point θ̇ = 0; L_z gives φ̇ there.
3. Cylinder formulas with ℓ = 2 m; T = L²/(2I) at fixed L.
4. `bang_bang_slew`; the body angle about axis 3 is 2 atan2(q_z, q_w) when the motion stays about that axis; for (c) the tilt of axis 3 is the angle between R e₃ and the inertial z axis; for (d) h_max = I ω_max and τ = 4IΘ/T².
5. `steady_precession_rate`, `fast_top_precession`, `simulate_rigid_body` with `torque=gravity`; (14) for (b); `nutation_turning_tilt` and `find_peaks` for (c).
6. Angular momentum accumulates as τt; I changes for the skater; rotor reaction torques; quaternion normalisation and frames; the energy and L_z invariants; energy dissipation at fixed L.
