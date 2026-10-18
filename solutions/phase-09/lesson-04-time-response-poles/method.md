# P09.04  Method

The route to each answer, without the final numbers. `worked-solutions.md` has the numbers for s = 15.

## 1. Reading the poles

Poles −ζω_n ± jω_n√(1 − ζ²). `overshoot_formula`, `peak_time_formula`, `settling_time_envelope`, 4/(ζω_n), `rise_time_exact`. Then `step_metrics(t, second_order_step(t, zeta, wn), 1.0)` on t = linspace(0, 6, 60001) and tabulate. The settling time disagrees most: the formulas bound the envelope.

## 2. Derivations

(a) A = 1; matching coefficients gives B = −1, C = −2σ; rewrite (s + 2σ)/((s + σ)² + ω_d²) as (s + σ)/(…) + (σ/ω_d)ω_d/(…); invert; combine cos and sin into one sine with phase φ, cos φ = ζ. (b) ẏ = (ω_n/√(1 − ζ²))e^{−σt} sin ω_d t; first zero at t > 0 is π/ω_d; y(t_p) = 1 + e^{−σπ/ω_d}; σ/ω_d = ζ/√(1 − ζ²). Inversion: square ln M = −πζ/√(1 − ζ²) and solve for ζ². (c) |y − 1| ≤ e^{−σt}/√(1 − ζ²); set equal to 0.02 and solve for t; once the envelope is inside the band the response is too, but not conversely. (d) T = (K/τ)/(s² + s/τ + K/τ): ω_n² = K/τ, 2ζω_n = 1/τ ⇒ ζω_n = 1/(2τ); underdamped when 1/τ² < 4K/τ. (e) Y_z = Y − (s/z)Y and sY ↔ ẏ (y(0) = 0); ẏ_z(0⁺) = ẏ(0⁺) − ÿ(0⁺)/z with ẏ(0⁺) = 0 and ÿ(0⁺) = ω_n².

## 3. Gimbal actuator

(a) Sector: angle from the negative real axis ≤ arccos ζ_min with ζ_min from (1) at 15 %; band: |Im s| ≥ π/0.06. (b) K = 1/(4ζ_min²τ) from `servo_gain_for_overshoot`; ω_d = ω_n√(1 − ζ²); t_p = π/ω_d. Any larger K increases overshoot, any smaller K increases t_p (ω_d falls): if t_p > 60 ms at the overshoot limit, no K works. (c) At ζ = ζ_min, ω_n = 1/(2ζ_min τ) and t_p = 2πζ_min τ/√(1 − ζ_min²) ≤ 0.06 gives the largest τ. (d) `simulate_step(servo_closed_loop(K, tau), t)` and `step_metrics`.

## 4. Filtered potentiometer

(a) T = K(τ_f s + 1)/(s(τs + 1)(τ_f s + 1) + K): the filter pole is a closed-loop zero at −1/τ_f. (b) `tf(num, den)` with `np.polymul`/`np.polyadd`; `np.roots`; `dominant_poles`; `residue_step_response` and `step_metrics` on linspace(0, 1.5, 6001). (c) Compare the new pair with the old one (−10 ± 13.65j at K = 14.31 without the filter): the filter has moved the pair, which is what changes the overshoot; the remote real pole and zero are almost a cancelling pair. (d) Bisection on K for overshoot = 10 %, then `step_metrics`.

## 5. Choosing the gain

(a) Loop τ ∈ {0.04, 0.05, 0.06} with K = 8.2; metrics on a fine grid. (b) Same with K = 14.31. (c) `np.arange(0, 1.5 + dt/2, dt)` with dt = 0.02 and 5 × 10⁻⁴; compare rise, peak time, overshoot and settling. Interpolation recovers rise and settling well; the peak time is quantised to the grid. (d) Weigh the nominal settling against its sensitivity to τ and the overshoot margin; a specification should state a robust figure with the conditions.

## 6. Right-half-plane zero

(a) From 2(e): ẏ_z(0⁺) = −ω_n²/z < 0 for z > 0, and y_z ends at +1, so it must first go negative. (b) `zero_step_closed_form(t, 0.5, 1.0, z)` on a fine grid; `step_metrics`; `t[np.argmin(y)]`; slope (y[1] − y[0])/dt. (c) T = n_G/(d_G + n_G): the zeros of T are those of n_G, independent of the gain (unless cancelled). (d) Elevator up: the tailplane lift is downward, the aircraft's total lift drops before the pitch rotation raises the angle of attack; altitude first falls. An aggressive altitude loop reacts to the initial sink with more elevator and can go unstable; the RHP zero limits the achievable bandwidth.
