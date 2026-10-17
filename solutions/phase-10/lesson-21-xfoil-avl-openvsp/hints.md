# P10.21  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `REQUIRED_SETTINGS["xfoil"]` lists what the reader insists on. XFOIL's polar file header carries the version, the airfoil name, xtrf, Mach, Re and Ncrit; look at what is missing from that list.
2. (a) α_L0 is a ratio of two correlated estimates: propagate with both variances and the covariance. (c) If cl₂ = (1 + ε) cl₁ and cl₁ = m α, the difference is ε m α.
3. Ladson gives both a free-transition and a grit-strip value at 6 × 10⁶. A cd below both means more laminar flow than the tunnel model had: what in XFOIL controls how long the flow stays laminar?
4. 15 mod 4 = 3. The worst cd ratio usually appears where transition moves fastest with α on the upper surface.
5. Run `vortex_panel` and `bl_drag` in loops; compare the spreads. A wraparound grit strip adds roughness drag of its own, not only an early transition.
6. The tunnel measured the physics; both codes are models of it. One-way coupling never feeds the boundary layer back into the inviscid solution.
