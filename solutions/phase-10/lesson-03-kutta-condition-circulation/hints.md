# P10.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Stagnation points satisfy sin(θ − α) = −Γ/(4πaV). The speed is largest where the circulation's velocity adds to the stream's, at the top for clockwise Γ.
2. (a) Multiply dw/dζ by ie^{iθ} on the circle and take the imaginary part with a minus sign. (b) w_z = w_ζ/(dz/dζ): a zero in the denominator needs a zero in the numerator. (c) c_l = ρVΓ/(½ρV²c). (d) The cross term between 2 sin θ and the constant is the only one that survives ∮ p sin θ dθ.
3. Γ = (surface speed) × (circumference) = 2πaU_s. If g > 1 the stagnation point is no longer on the surface.
4. `fit_karman_trefftz`, `section_shape`, `circle_geometry`, then c_l = 8πa sin(α_x + β)/c with α measured from the x-axis; convert to the chord line with the inclination φ.
5. Joukowski's dz/dζ has a simple zero at ζ = b; the Karman–Trefftz derivative behaves like r^{n−1} with a non-integer power.
6. A Kutta condition at the trailing edge assumes the boundary layer is thin there. Is it, at Re = 2 × 10⁶?
