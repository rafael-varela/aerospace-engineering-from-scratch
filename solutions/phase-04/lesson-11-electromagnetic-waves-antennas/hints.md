# P04.11  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. EIRP = P_t G_t; then (3) with EIRP in place of P; then A_e = λ²/(4π) for G = 1. Convert to dBm with 10 log₁₀(P/1 mW).
2. (a) Curl of Faraday's law and the identity ∇ × ∇ × E = ∇(∇·E) − ∇²E. (b) ∫₀^π sin³θ dθ = 4/3 and ∫ sin⁵θ dθ = 16/15. (c) The transmitter's intensity at the receiver is G_t P_t/(4πd²); the receiver collects it over A_e,r. (d) The law of cosines with the angle at the station equal to 90° + e.
3. Both the light intensity and gravity fall as 1/r². For the transverse component, maximise cos²α sin α.
4. `dish_gain(D, C_LIGHT/2.25e9)`, `effective_aperture`; `Link(frequency=2.25e9, P_t=2.0, G_t_dBi=6.0, G_r_dBi=..., T_sys=300.0, other_losses_dB=4.0, bit_rate=1e6)`. The bit rate scales with 10^{(margin − 3)/10}.
5. A lambda `f = lambda t: dipole_pattern(t, L)` passes your length to the integrators. The feed sees the current at the centre, I₀ sin(kL/2).
6. Put a number on every claim: gain in dBi, beamwidth in degrees, cost per decibel.
