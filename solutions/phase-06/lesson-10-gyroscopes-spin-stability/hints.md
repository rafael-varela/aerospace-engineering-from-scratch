# P06.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `nutation(800, 1200, w3, 0.02)`; the ground sees the axis move in space, not in the body.
2. Write H's components along and across the symmetry axis; for (c), add the two equations of (8) and dot each with its own rate; for (e), a stiff coupling makes the sphere part of the body.
3. τ = F × offset; H = I₃ω₃; the spin axis also nutates around H; from rest, angle = ½(τ/I)t².
4. `gyroscope_torque(I, tau)` with a starting quaternion that tilts the body about inertial x; fit the unwrapped azimuth of H_i.
5. `kane_nutation_mode` in a loop (or `scipy.optimize.minimize_scalar` on log c); `simulate_with_damper`, `fit_growth_rate`, `spin_axis_inertial`, `damper_energy`; ln(0.2/10⁻³)/σ.
6. Rigid versus dissipative; the major-axis rule; burn time against e-folding time; the sign of the energy-sink rate; (7) against (9); illustrative inertia and damper.
