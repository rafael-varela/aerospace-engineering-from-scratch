# P03.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything is a function of x = (hω)². The determinant check is (1 − x/2)² + x(1 − x/4). The stability limit hω = 2 corresponds to π steps per period.
2. (a) A kick has Jacobian [[I, 0], [hA, I]] with A = ∂a/∂q; compute MᵀJM and watch the term hA − hAᵀ appear. (b) q × (q f) = 0 and p × p = 0. (c) Write the three sub-maps of Verlet as K_{h/2} D_h K_{h/2} and reverse the order with −h; for symplectic Euler, K_h D_h followed by K_h D_h is K_h (D_h K_h D_h) K_h... regroup.
3. Angular rate = speed/radius; v_p/v_a = r_a/r_p by angular momentum. The waste fraction is 1 − (apogee count)/(perigee count). For (c) remember that Verlet's swing scales as (hω_p)² and RK4's error as (hω_p)⁴: at hω_p = 0.3 the fourth power is much smaller.
4. `perigee_state(6300e3 / (1 - e), e)`; `energy_drift` with 50 orbits and 200 steps; `drift_slope`; `kepler_exact(50 * P, q0, p0)`. The perigee passage at higher e is shorter in *time* for the same hω_p, so fewer steps resolve it.
5. The crossover is where (slope × N) = swing. Yoshida's cost is three Verlet steps per step. Position error and energy error measure different things: phase along the orbit against size of the orbit.
6. Ask two questions: is the force conservative, and is the propagation long enough for a secular drift to matter? Burns, drag and events answer the first; the orbit count answers the second.
