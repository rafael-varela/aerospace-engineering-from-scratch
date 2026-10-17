# P08.04  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The manometer gives Δp = (ρ_w − ρ_air) g Δh, and the pitot-static difference is ½ρV². For the uncertainty, V ∝ Δh^{1/2}: take logarithms and differentiate.
2. (a) Along a streamline the acceleration of a particle is V ∂V/∂s in steady flow. (b) Expand (1 + x)^{3.5} with x = 0.2M² to second order. (c) EAS keeps the dynamic pressure, and the wing's lift and stall depend on dynamic pressure.
3. ρgQ times the head is the power the liquid receives; the head is dominated by whichever of pressure, velocity and height changes most per unit weight.
4. The airspeed chain needs the pressure, temperature and density at the pressure altitude; the port error changes the sensed impact pressure by −c_p q. For (c), the IAS is a monotonic function of TAS, so a bracketing root finder works.
5. The ratio q_c/q depends only on M (write it out); so does the Bernoulli error. For the contraction, continuity gives V_in = V_test/CR.
6. Each condition of Bernoulli has a place where it fails on an aeroplane: work is done (propeller), viscosity matters (boundary layer), the flow separates and mixes (wake), the density changes (high Mach number), the flow changes in time (gust).
