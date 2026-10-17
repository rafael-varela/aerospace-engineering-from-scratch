# P04.12  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. A ray entering the end face at θ_a refracts to θ_core inside; it then meets the wall at 90° − θ_core, which must be at least θ_c. Use cos θ_c = √(1 − (n₂/n₁)²).
2. (a) A star at angle α images at f tan α. (b) Write each recorded position as true position + photon scatter + pixel rounding, and assume the true positions are spread so that the rounding is uniform over a pixel. (c) Each side of the bridge is a divider; subtract the two divider outputs. (d) Expand (1 + ε)^k = 1 + kε + k(k − 1)ε²/2 with ε = (γ − 1)M²/2, and use ρV² = γpM².
3. V = V_s R_T/(R_top + R_T); differentiate with respect to R_T and then multiply by dR_T/dT; maximise over R_top. Self-heating error = power/dissipation constant.
4. `StarTracker(f=..., aperture=f/1.4, n_pixels=2048, pitch=3.45e-6, m_limit=5.5, t_exp=0.05)`; its properties give the scale, FOV and star count. `simulate_centroid_errors(Ne, 1.0, 3.0, 2000, size=..., seed=s)`.
5. `quarter_bridge(3.3, 1500e-6)`, `full_bridge(3.3, 1500e-6)`; the gain is 1.6 V divided by the bridge output; the strain per count is the strain that moves the ADC input by 3.3 V/4096.
6. Every verdict should point to a number already in the lesson: the pixel-locking bias, the window-and-read-noise table, the 1/√n of the attitude error, and what one star can and cannot tell you about attitude.
