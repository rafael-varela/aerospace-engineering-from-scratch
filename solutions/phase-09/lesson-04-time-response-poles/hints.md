# P09.04  Hints

One hint per problem. Read only the one you need, then `method.md` if the hint is not enough.

1. σ = ζω_n and ω_d = ω_n√(1 − ζ²). The envelope bound and 4/(ζω_n) bound the envelope, not the response.

2. (a) Write the quadratic as (s + σ)² + ω_d² and split G/s into A/s + (Bs + C)/((s + σ)² + ω_d²). (b) Differentiate y and find the first positive zero of ẏ. (c) |sin| ≤ 1. (d) Read ζ and ω_n from τs² + s + K after dividing by τ. (e) G_z = G − (s/z)G; the initial value theorem gives ẏ(0⁺) = lim s·sY.

3. Overshoot fixes a minimum ζ (a sector around the negative real axis); peak time fixes a minimum ω_d (outside a horizontal band). With proportional control the real part is fixed at −1/(2τ).

4. The filter pole becomes a closed-loop zero, and the characteristic polynomial is s(τs + 1)(τ_f s + 1) + K. The dominance ratio says nothing about whether the pair itself has moved.

5. Use `servo_closed_loop`, `second_order_params` and `second_order_step` on two grids. A settling time that depends on the first peak just touching the band is fragile.

6. ẏ(0) = 0 for the standard system, so the slope comes from the −(1/z)ÿ term at t = 0⁺. In T = n_G/(d_G + n_G), the numerator is n_G whatever the gain.
