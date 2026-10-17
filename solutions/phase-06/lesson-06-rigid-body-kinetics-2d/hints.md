# P06.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Assume it rolls, find the friction that needs, and compare with μ_s N; (4) gives both (b) and the steepest ramp in (c).
2. Differentiate v_s = V − ωR and substitute both equations of (5); for the energy, F dt = −M_eff dv_s turns a time integral into an integral over v_s.
3. The canister is its own "axle": M = m and N = mg. Backspin makes ω₀ negative, so v_s0 = V₀ − ω₀R is larger than V₀; the sign of V_f decides whether it returns.
4. Check J against ½ μ_k N_max t_r before choosing a branch of (8); `spin_up_closed_form` and `simulate_touchdown` with `ramp_load(N_max, 0.08)`.
5. `max_brake_torque`; a brake function `lambda t: T if t >= 0.3 else 0.0`; read `res.transitions`, `res.E_brake` and `res.E_friction`; d = V₀²/(2a).
6. Energy per unit of tread area; (1 − f)² and the sign of v_s0; the moment equation about an accelerating point; terminal events; what a sudden load release does to a spring; which results need μ and N(t), and which need only momentum.
