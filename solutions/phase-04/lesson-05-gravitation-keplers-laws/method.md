# P04.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. r = R + (500 + 10s) km. 2. g = μ/r², v = √(μ/r), T = 2πr/v, ε = −μ/2r. 3. v_esc = √2 v; Δv = (√2 − 1)v. 4. Check v²/r = g and v²/2 = μ/2r.

## 2.
(a) h = r_p v_p = r_a v_a; substitute v_a into ε at both apsides; solve for v_p²; put back into ε; r_p + r_a = 2a; then ε = v²/2 − μ/r everywhere. (b) h² = r_p² v_p² with v_p² from (a) and r_p = a(1 − e), r_a = a(1 + e); solve for e using ε = −μ/2a. e = 1 is ε = 0. (c) dA/dt = h/2 constant ⇒ πab = (h/2)T; substitute b and h. (d) Kick: Δ(r × v) = r × (h/2) f r = 0. Drift: (r + hv) × v = r × v. RK4 combines slopes evaluated at different intermediate states; the update of v is not parallel to the final r, so the cross product changes at O(h⁵).

## 3.
(a) r = (μT²/4π²)^{1/3}; altitude r − R; v = √(μ/r). (b) g = μ/R², v_esc = √(2μ/R). (c) T = 2π√(a³/μ); compare with 88 642.66 s; an orbit faster than the surface turns overtakes it, so Phobos moves west-to-east across the sky faster than the stars move east-to-west.

## 4.
(a) `perigee_state(a, e)`; altitudes a(1 ∓ e) − R; `vis_viva` at both. (b) Loop over steps per orbit: `propagate_leapfrog(a, e, 10, spo)` for energy and h, `propagate_leapfrog(a, e, 1, spo)` interpolated onto `propagate_solve_ivp(a, e, 1)` times for position; `convergence_order`. (c) N = 1000 √(ΔE₁₀₀₀/10⁻⁶); round up and rerun. (d) `period_from_trajectory` on a five-orbit run; subtract `orbital_period(a)`.

## 5.
(a) `hohmann_transfer(R + (300 + s) km, 26 560 km)`; e_t = (r₂ − r₁)/(r₂ + r₁). (b) v = v_c(r₁) + Δv₁ along +y at (r₁, 0); `orbital_elements` returns r_a. (c) Plane-change cost ∝ the speed at which it is done; at apoapsis of the transfer the speed is lowest, so combine the plane change with Δv₂ (the vector sum is less than the two burns done separately).

## 6.
For each case: list the missing physics (J₂, drag, third bodies, solar pressure), estimate its size over the stated time, compare with the required accuracy; for numerics, compare the leapfrog's phase error and its step efficiency at high eccentricity.
